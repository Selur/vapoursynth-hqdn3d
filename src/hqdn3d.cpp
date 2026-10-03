/*
    HQDN3D 1.10 for Vapoursynth

    Copyright (C) 2003 Daniel Moreno <comac@comac.darktech.org>
    Avisynth port (C) 2005 Loren Merritt <lorenm@u.washington.edu>
    Vapoursynth port (C) 2017 Martin Güthle  <mguethle@xunit.de>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <vector>

#include <VapourSynth4.h>
#include <VSHelper4.h>

typedef struct Hqdn3dData {
    VSNode *clip;
    const VSVideoInfo *vi;

    double lumSpac;
    double chromSpac;
    double lumTmp;
    double chromTmp;
    int restartLap;

    // Pixels are processed as 24 bit fixed point numbers: an 8 bit sample
    // with 16 fractional bits. Samples of other bit depths are shifted into
    // that scale on load and back on store, so the filter core and the
    // strength parameters are independent of the bit depth of the clip.
    int inShift;   // 24 - bitsPerSample
    int peak;      // (1 << bitsPerSample) - 1
    // The coefficient tables are indexed by the pixel difference at a
    // resolution of 1 / (1 << lutBits) of an 8 bit step. 4 keeps the
    // original 8 bit tables; from 12 bit on one step equals one source LSB.
    int lutBits;
    int lutShift;  // 16 - lutBits
    std::vector<int> coefs[4];
    uint16_t *prevFrame[3] = { nullptr, nullptr, nullptr };
    unsigned int *prevLine[3] = { nullptr, nullptr, nullptr };
    bool process[3];
    int last_frame = -1; // the last frame returned by hqdn3d
} Hqdn3dData;

static void
VS_CC hqdn3dFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    (void)core;

    Hqdn3dData *d = (Hqdn3dData *)instanceData;

    for (int p = 0; p < d->vi->format.numPlanes; p++) {
        free(d->prevLine[p]);
        free(d->prevFrame[p]);
    }

    vsapi->freeNode(d->clip);
    delete d;
}

static inline unsigned int
LowPassMul(unsigned int pMul, unsigned int cMul, const int* coef, const int lutShift) {
    // The 1 << 24 bias keeps the difference positive before the shift, the
    // remainder rounds it to the nearest table entry. For lutShift == 12
    // this is the original 0x10007FF.
    const unsigned int roundConvolution = (1u << 24) + (1u << (lutShift - 1)) - 1;
    int d = (static_cast<int>(pMul - cMul) + roundConvolution) >> lutShift;
    return cMul + coef[d];
}

static inline unsigned int
loadPixel(const uint8_t sample, const int inShift, const int peak) {
    (void)peak;
    return static_cast<unsigned int>(sample) << inShift;
}

static inline unsigned int
loadPixel(const uint16_t sample, const int inShift, const int peak) {
    // Samples above the nominal range would index the coefficient tables out
    // of bounds, so clamp them.
    return static_cast<unsigned int>(std::min<int>(sample, peak)) << inShift;
}

template <typename PixelT>
static void
deNoise(
      const uint8_t * srcPlane
    , uint16_t * prevPlane
    , unsigned int * prevLine
    , uint8_t * tarPlane
    , const int frameWidth
    , const int frameHeight
    , const int srcStride
    , const int tarStride
    , const int *coefsHorizontal
    , const int *coefsVertical
    , const int *coefsTemporal
    , const bool isFirstFrame
    , const int inShift
    , const int peak
    , const int lutShift
) {
    static const unsigned int ROUND_LINE  = 0x1000007F;
    static const unsigned int SHIFT_LINE  = 8;
    // Rounds like the original 0x10007FFF: ties go down.
    const unsigned int roundPixel = (1u << (inShift - 1)) - 1;

    const int srcPitch = srcStride / static_cast<int>(sizeof(PixelT));
    const int tarPitch = tarStride / static_cast<int>(sizeof(PixelT));

    for (int row = 0; row < frameHeight; ++row) {
        const PixelT *srcRow = reinterpret_cast<const PixelT *>(srcPlane) + row * srcPitch;
        PixelT *tarRow = tarPlane
            ? reinterpret_cast<PixelT *>(tarPlane) + row * tarPitch
            : nullptr;
        /* gcc assume prevPixel might be used in an uninitialized way, but
         * it's not. So feel free to use any other value
         */
        unsigned int prevPixel = 0;
        for (int col = 0; col < frameWidth; ++col) {
            const unsigned int curPixel = loadPixel(srcRow[col], inShift, peak);
            // Correlate current pixel with previous pixel
            prevPixel = col == 0
                ? curPixel
                : LowPassMul(
                      prevPixel
                    , curPixel
                    , coefsHorizontal
                    , lutShift
                );
            // Correlate previous line with previous pixel
            prevLine[col] = row == 0
                ? prevPixel
                : LowPassMul(
                      prevLine[col]
                    , prevPixel
                    , coefsVertical
                    , lutShift
                );
            unsigned int resPix;
            if (isFirstFrame) {
                resPix = prevLine[col];
            } else {
                // Correlate vertical result with previous result frame pixel
                resPix = LowPassMul(
                      static_cast<unsigned int>(prevPlane[row * frameWidth + col]) << SHIFT_LINE
                    , prevLine[col]
                    , coefsTemporal
                    , lutShift
                );
            }

            // The temporal state keeps 16 bits: an 8 bit sample with 8
            // fractional bits, independent of the clip's bit depth.
            prevPlane[row * frameWidth + col]
                = static_cast<uint16_t>(((resPix + ROUND_LINE) >> SHIFT_LINE) & 0xFFFF);
            if (tarRow)
                tarRow[col] = static_cast<PixelT>((resPix + roundPixel) >> inShift);
        }
    }
}

static void filterFrame(
      const VSFrame *srcFrame
    , VSFrame *newFrame
    , const bool isFirstFrame
    , Hqdn3dData *usrData
    , const VSAPI *vsapi
) {
    const VSVideoFormat *srcFrameFmt = vsapi->getVideoFrameFormat(srcFrame);

    for (int plane = 0; plane < srcFrameFmt->numPlanes; plane++) {
        if (!usrData->process[plane])
            continue;

        auto deNoiseFn = srcFrameFmt->bytesPerSample == 1
            ? deNoise<uint8_t>
            : deNoise<uint16_t>;

        deNoiseFn(
              vsapi->getReadPtr(srcFrame, plane)
            , usrData->prevFrame[plane]
            , usrData->prevLine[plane]
            , newFrame ? vsapi->getWritePtr(newFrame, plane) : nullptr
            , vsapi->getFrameWidth(srcFrame, plane)
            , vsapi->getFrameHeight(srcFrame, plane)
            , static_cast<int>(vsapi->getStride(srcFrame, plane))
            , newFrame ? static_cast<int>(vsapi->getStride(newFrame, plane)) : 0
            , usrData->coefs[plane == 0 ? 0 : 2].data() // Y or U/V
            , usrData->coefs[plane == 0 ? 0 : 2].data() // Y or U/V
            , usrData->coefs[plane == 0 ? 1 : 3].data() // Y or U/V
            , isFirstFrame
            , usrData->inShift
            , usrData->peak
            , usrData->lutShift
        );
    }

}

static const VSFrame *VS_CC hqdn3dGetFrame(
      int n
    , int activationReason
    , void *instanceData
    , void **frameData
    , VSFrameContext *frameCtx
    , VSCore *core
    , const VSAPI *vsapi
) {
    (void)frameData;

    // Get the user data
    Hqdn3dData * usrData = reinterpret_cast<Hqdn3dData *>(instanceData);

    if (activationReason == arInitial) {
        // if we skip some frames, filter the gap anyway
        if (n > usrData->last_frame + 1 &&
            n - usrData->last_frame <= usrData->restartLap + 1 &&
            usrData->last_frame >= 0) {

            for (int i = usrData->last_frame + 1; i < n; i++) {
                vsapi->requestFrameFilter(i, usrData->clip, frameCtx);
            }
        // if processing out of sequence, filter several previous frames to minimize seeking problems
        } else if (n != usrData->last_frame + 1) {
            int sn = std::max(0, n - usrData->restartLap);

            for (int i = sn + 1; i < n; i++)
                vsapi->requestFrameFilter(i, usrData->clip, frameCtx);
        }

        vsapi->requestFrameFilter(n, usrData->clip, frameCtx);

        return nullptr;
    }
    if (activationReason != arAllFramesReady) {
        return nullptr;
    }

    // if we skip some frames, filter the gap anyway
    if (n > usrData->last_frame + 1 &&
        n - usrData->last_frame <= usrData->restartLap + 1 &&
        usrData->last_frame >= 0) {

        for (int i = usrData->last_frame + 1; i < n; i++) {
            const VSFrame *f = vsapi->getFrameFilter(i, usrData->clip, frameCtx);

            filterFrame(f, nullptr, false, usrData, vsapi);

            vsapi->freeFrame(f);
        }
    // if processing out of sequence, filter several previous frames to minimize seeking problems
    } else if (n != usrData->last_frame + 1) {
        int sn = std::max(0, n - usrData->restartLap);

        for (int i = sn + 1; i < n; i++) {
            const VSFrame *f = vsapi->getFrameFilter(i, usrData->clip, frameCtx);

            filterFrame(f, nullptr, i == sn + 1, usrData, vsapi);

            vsapi->freeFrame(f);
        }
    }


    // Get current frame
    const VSFrame * srcFrame =
        vsapi->getFrameFilter(n, usrData->clip, frameCtx);


    // Create target frame
    const VSFrame *plane_src[3] = {
        usrData->process[0] ? nullptr : srcFrame,
        usrData->process[1] ? nullptr : srcFrame,
        usrData->process[2] ? nullptr : srcFrame
    };
    int planes[3] = { 0, 1, 2 };

    VSFrame *newFrame = vsapi->newVideoFrame2(
          &usrData->vi->format
        , usrData->vi->width
        , usrData->vi->height
        , plane_src
        , planes
        , srcFrame
        , core
    );

    filterFrame(srcFrame, newFrame, n == 0, usrData, vsapi);

    vsapi->freeFrame(srcFrame);

    usrData->last_frame = n;

    return newFrame;
}


// Create the plugin
static void VS_CC hqdn3dCreate(
      const VSMap *in
    , VSMap *out
    , void *userData
    , VSCore *core
    , const VSAPI *vsapi
) {
    (void)userData;

    Hqdn3dData d;
    int err;

    d.lumSpac    = vsapi->mapGetFloat(in, "lum_spac",    0, &err);
    if (err) {
        d.lumSpac = 4.0;
    } else if (d.lumSpac < 0 || d.lumSpac > 255) {
        vsapi->mapSetError(out, "Hqdn3d: lum_spac must be between 0 and 255 (inclusive).");
        return;
    }

    d.chromSpac  = vsapi->mapGetFloat(in, "chrom_spac",  0, &err);
    if (err) {
        d.chromSpac = .75 * d.lumSpac;
    } else if (d.chromSpac < 0 || d.chromSpac > 255) {
        vsapi->mapSetError(out, "Hqdn3d: chrom_spac must be between 0 and 255 (inclusive).");
        return;
    }

    d.lumTmp     = vsapi->mapGetFloat(in, "lum_tmp",     0, &err);
    if (err) {
        d.lumTmp = 1.5 * d.lumSpac;
    } else if (d.lumTmp < 0 || d.lumTmp > 255) {
        vsapi->mapSetError(out, "Hqdn3d: lum_tmp must be between 0 and 255 (inclusive).");
        return;
    }

    d.chromTmp   = vsapi->mapGetFloat(in, "chrom_tmp",   0, &err);
    if (err) {
        d.chromTmp = (d.lumSpac == 0) ? d.chromSpac * 1.5
                                      : d.lumTmp * d.chromSpac / d.lumSpac;
    } else if (d.chromTmp < 0 || d.chromTmp > 255) {
        vsapi->mapSetError(out, "Hqdn3d: chrom_tmp must be between 0 and 255 (inclusive).");
        return;
    }

    d.restartLap = vsh::int64ToIntS(vsapi->mapGetInt(in, "restart_lap", 0, &err));
    if (err)
        d.restartLap = std::max(2
            , static_cast<int>(1 + std::max(d.lumTmp, d.chromTmp)));


    d.clip = vsapi->mapGetNode(in, "clip", 0, nullptr);
    d.vi = vsapi->getVideoInfo(d.clip);

    if (!vsh::isConstantVideoFormat(d.vi) ||
        d.vi->format.colorFamily == cfRGB ||
        d.vi->format.sampleType != stInteger ||
        d.vi->format.bitsPerSample < 8 ||
        d.vi->format.bitsPerSample > 16) {

        vsapi->mapSetError(out, "Hqdn3d: input clip must be 8-16 bit integer, not RGB, and it must have constant format and dimensions.");
        vsapi->freeNode(d.clip);
        return;
    }

    const int bitsPerSample = d.vi->format.bitsPerSample;
    d.inShift  = 24 - bitsPerSample;
    d.peak     = (1 << bitsPerSample) - 1;
    d.lutBits  = std::max(4, bitsPerSample - 8);
    d.lutShift = 16 - d.lutBits;


    d.lumSpac   = std::min(254.9, d.lumSpac);
    d.chromSpac = std::min(254.9, d.chromSpac);
    d.lumTmp    = std::min(254.9, d.lumTmp);
    d.chromTmp  = std::min(254.9, d.chromTmp);

    // Calculate the coefficients. The strengths are always on the 8 bit
    // scale, the tables map a pixel difference (in 1 / lutSteps of an 8 bit
    // step) to the correction in the 24 bit internal scale.
    const int lutSteps  = 1 << d.lutBits;
    const int lutCenter = 256 * lutSteps;
    for (auto const &cc : {
          std::make_pair(0, d.lumSpac)
        , std::make_pair(1, d.lumTmp)
        , std::make_pair(2, d.chromSpac)
        , std::make_pair(3, d.chromTmp)
    } ) {
        const double gamma = std::log(0.25) / std::log(1.0 - cc.second / 255.0 - 0.00001);
        std::vector<int> &table = d.coefs[cc.first];
        table.assign(512 * lutSteps, 0);
        for (int i = -255 * lutSteps; i < 256 * lutSteps; ++i) {
            const double simil = std::max(0.0, 1.0 - std::abs(i) / (lutSteps * 255.0));
            const double c = std::pow(simil, gamma) * 65536.0 * i / lutSteps;
            table[lutCenter + i]
                = static_cast<int>(c < 0 ? c - 0.5 : c + 0.5);
        }
    }

    // According to the documentation, 0 strength means no processing.
    d.process[0] = d.lumSpac != 0 || d.lumTmp != 0;
    d.process[1] = d.process[2] = d.chromSpac != 0 || d.chromTmp != 0;


    for (int p = 0; p < d.vi->format.numPlanes; p++) {
        int width = d.vi->width;
        int height = d.vi->height;

        if (p) {
            width >>= d.vi->format.subSamplingW;
            height >>= d.vi->format.subSamplingH;
        }

        d.prevFrame[p] = (uint16_t *)malloc(width * height * sizeof(uint16_t));
        d.prevLine[p] = (unsigned int *)malloc(width * sizeof(unsigned int));
    }


    Hqdn3dData *data = new Hqdn3dData(d);

    // The filter keeps temporal state, so frames must be requested one at a
    // time (fmFrameState), and it may request earlier frames (rpGeneral).
    VSFilterDependency deps[] = { { data->clip, rpGeneral } };
    vsapi->createVideoFilter(
          out
        , "Hqdn3d"
        , data->vi
        , hqdn3dGetFrame
        , hqdn3dFree
        , fmFrameState
        , deps
        , 1
        , data
        , core
    );
}

VS_EXTERNAL_API(void) VapourSynthPluginInit2(
      VSPlugin *plugin
    , const VSPLUGINAPI *vspapi
) {
    vspapi->configPlugin(
          "com.vapoursynth.hqdn3d"
        , "hqdn3d"
        , "HQDn3D port as used in avisynth/mplayer"
        , VS_MAKE_VERSION(1, 1)
        , VAPOURSYNTH_API_VERSION
        , 0
        , plugin
    );
    vspapi->registerFunction(
          "Hqdn3d"
        , "clip:vnode;"
          "lum_spac:float:opt;"
          "chrom_spac:float:opt;"
          "lum_tmp:float:opt;"
          "chrom_tmp:float:opt;"
          "restart_lap:int:opt;"
        , "clip:vnode;"
        , hqdn3dCreate
        , nullptr
        , plugin
    );
}
