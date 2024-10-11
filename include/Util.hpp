#pragma once
#include <AK/SoundEngine/Common/AkCommonDefs.h>
#include <AK/SoundEngine/Common/AkSimd.h>

namespace Wpe
{
    constexpr AkUInt32 DEFAULT_SIMD_WIDTH = 4;

    inline void OverlapAdd(AkSampleType* io_pOrg, const AkSampleType* in_pLayer, AkUInt32 in_uNumSamples)
    {
        AkUInt32 uVectorIterations = in_uNumSamples / DEFAULT_SIMD_WIDTH;
        while (uVectorIterations--)
        {
            AKSIMD_V4F32 vLayer = AKSIMD_LOAD_V4F32(in_pLayer);
            AKSIMD_V4F32 vOrg = AKSIMD_LOAD_V4F32(io_pOrg);
            AKSIMD_V4F32 vOut = AKSIMD_ADD_V4F32(vOrg, vLayer);
            AKSIMD_STORE_V4F32(io_pOrg, vOut);

            in_pLayer += DEFAULT_SIMD_WIDTH;
            io_pOrg += DEFAULT_SIMD_WIDTH;
        }

        const auto uRemainingSamples = in_uNumSamples % DEFAULT_SIMD_WIDTH;
        for (AkUInt32 i = 0; i < uRemainingSamples; ++i)
        {
            *io_pOrg += *in_pLayer;
            ++io_pOrg;
            ++in_pLayer;
        }
    }

    // https://signalsmith-audio.co.uk/writing/2021/cheap-energy-crossfade/#fractional-power-law
    inline void EnergyCrossfadePair(AkSampleType x, AkSampleType& fadeIn, AkSampleType& fadeOut)
    {
        AkSampleType x2 = 1 - x;
        AkSampleType A = x * x2;
        AkSampleType B = A * (1 + 1.4186f * A);
        AkSampleType C = (B + x), D = (B + x2);
        fadeIn = C * C;
        fadeOut = D * D;
    }

    inline void EnergyCrossfadePair(AKSIMD_V4F32 x, AKSIMD_V4F32& fadeIn, AKSIMD_V4F32& fadeOut)
    {
        AKSIMD_V4F32 x2 = AKSIMD_SUB_V4F32(AKSIMD_SET_V4F32(1.0f), x);
        AKSIMD_V4F32 A = AKSIMD_MUL_V4F32(x, x2);
        AKSIMD_V4F32 B = AKSIMD_MUL_V4F32(A, AKSIMD_ADD_V4F32(AKSIMD_SET_V4F32(1.0f), AKSIMD_MUL_V4F32(AKSIMD_SET_V4F32(1.4186f), A)));
        AKSIMD_V4F32 C = AKSIMD_ADD_V4F32(B, x);
        AKSIMD_V4F32 D = AKSIMD_ADD_V4F32(B, x2);
        fadeIn = AKSIMD_MUL_V4F32(C, C);
        fadeOut = AKSIMD_MUL_V4F32(D, D);
    }

    inline void Crossfade(AkSampleType* in_pFadeOut, AkSampleType* in_pFadeIn, AkSampleType* io_pDst, AkUInt32 in_uNumSamples)
    {
        AkUInt32 uVectorIterations = in_uNumSamples / DEFAULT_SIMD_WIDTH;
        for (AkUInt32 i = 0; i < uVectorIterations; ++i)
        {
            AKSIMD_V4F32 vFadeIn, vFadeOut;
            auto vNumSamples = AKSIMD_SET_V4F32(in_uNumSamples);
            auto fBaseIndex = DEFAULT_SIMD_WIDTH * i;
            auto vIndex = AKSIMD_SETV_V4F32(fBaseIndex+3, fBaseIndex+2, fBaseIndex+1, fBaseIndex);
            EnergyCrossfadePair(AKSIMD_DIV_V4F32(vIndex, vNumSamples), vFadeIn, vFadeOut);
            auto vDst = AKSIMD_ADD_V4F32(
                AKSIMD_MUL_V4F32(AKSIMD_LOAD_V4F32(in_pFadeIn), vFadeIn),
                AKSIMD_MUL_V4F32(AKSIMD_LOAD_V4F32(in_pFadeOut), vFadeOut)
            );
            AKSIMD_STORE_V4F32(io_pDst, vDst);
            in_pFadeOut += DEFAULT_SIMD_WIDTH;
            in_pFadeIn += DEFAULT_SIMD_WIDTH;
            io_pDst += DEFAULT_SIMD_WIDTH;
        }

        auto uRemainingSamples = in_uNumSamples % DEFAULT_SIMD_WIDTH;
        auto uOffset = DEFAULT_SIMD_WIDTH * uVectorIterations;
        for (AkUInt32 i = 0; i < uRemainingSamples; ++i)
        {
            AkSampleType fFadeIn, fFadeOut;
            auto uIndex = uOffset + i;
            EnergyCrossfadePair(static_cast<AkSampleType>(uIndex) / in_uNumSamples, fFadeIn, fFadeOut);
            io_pDst[i] = (fFadeIn * in_pFadeIn[uIndex]) + (fFadeOut * in_pFadeOut[uIndex]);
        }
    }
}
