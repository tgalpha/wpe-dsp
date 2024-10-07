#pragma once
#include <AK/SoundEngine/Common/AkSimd.h>
#include <AK/SoundEngine/Common/AkCommonDefs.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

namespace Wpe
{
    enum FilterTypes
    {
        LowPass,
        HighPass,
        BandPass,
        Notch,
        Peak,
        LowShelf,
        HighShelf,
        AllPass,
    };

    inline void ComputeBiquadCoefficients(AkUInt32 in_uSampleRate,
                                          FilterTypes in_eFilterType,
                                          AkReal32 in_fFreq,
                                          AkReal32 in_fQ,
                                          AkReal32 in_fGainDB,
                                          float* out_pB0, float* out_pB1, float* out_pB2,
                                          float* out_pA1, float* out_pA2)
    {
        const auto fOmega = 2.0f * M_PI * in_fFreq / in_uSampleRate;
        const auto fSin = sinf(fOmega);
        const auto fCos = cosf(fOmega);
        const auto fAlpha = fSin / (2.0f * in_fQ);
        const auto fGainAbs = powf(10.f, in_fGainDB * 0.025f);
        const auto fBeta = sqrtf(fGainAbs * 2);
        auto fA0 = 0.0f;

        switch (in_eFilterType)
        {
        case LowPass:
            *out_pB0 = (1.0f - fCos) / 2.0f;
            *out_pB1 = 1.0f - fCos;
            *out_pB2 = (1.0f - fCos) / 2.0f;
            fA0 = 1.0f + fAlpha;
            *out_pA1 = -2.0f * fCos;
            *out_pA2 = 1.0f - fAlpha;
            break;
        case HighPass:
            *out_pB0 = (1.0f + fCos) / 2.0f;
            *out_pB1 = -(1.0f + fCos);
            *out_pB2 = (1.0f + fCos) / 2.0f;
            fA0 = 1.0f + fAlpha;
            *out_pA1 = -2.0f * fCos;
            *out_pA2 = 1.0f - fAlpha;
            break;
        case BandPass:
            *out_pB0 = fAlpha;
            *out_pB1 = 0.0f;
            *out_pB2 = -fAlpha;
            fA0 = 1.0f + fAlpha;
            *out_pA1 = -2.0f * fCos;
            *out_pA2 = 1.0f - fAlpha;
            break;
        case Peak:
            *out_pB0 = 1.0f + fAlpha * fGainAbs;
            *out_pB1 = -2.0f * fCos;
            *out_pB2 = 1.0f - fAlpha * fGainAbs;
            fA0 = 1.0f + fAlpha / fGainAbs;
            *out_pA1 = -2.0f * fCos;
            *out_pA2 = 1.0f - fAlpha / fGainAbs;
            break;
        case Notch:
            *out_pB0 = 1.0f;
            *out_pB1 = -2.0f * fCos;
            *out_pB2 = 1.0f;
            fA0 = 1.0f + fAlpha;
            *out_pA1 = -2.0f * fCos;
            *out_pA2 = 1.0f - fAlpha;
            break;
        case LowShelf:
            *out_pB0 = fGainAbs * ((fGainAbs + 1.0f) - (fGainAbs - 1.0f) * fCos + fBeta * fSin);
            *out_pB1 = 2.0f * fGainAbs * ((fGainAbs - 1.0f) - (fGainAbs + 1.0f) * fCos);
            *out_pB2 = fGainAbs * ((fGainAbs + 1.0f) - (fGainAbs - 1.0f) * fCos - fBeta * fSin);
            fA0 = (fGainAbs + 1.0f) + (fGainAbs - 1.0f) * fCos + fBeta * fSin;
            *out_pA1 = -2.0f * ((fGainAbs - 1.0f) + (fGainAbs + 1.0f) * fCos);
            *out_pA2 = (fGainAbs + 1.0f) + (fGainAbs - 1.0f) * fCos - fBeta * fSin;
            break;
        case HighShelf:
            *out_pB0 = fGainAbs * ((fGainAbs + 1.0f) + (fGainAbs - 1.0f) * fCos + fBeta * fSin);
            *out_pB1 = -2.0f * fGainAbs * ((fGainAbs - 1.0f) + (fGainAbs + 1.0f) * fCos);
            *out_pB2 = fGainAbs * ((fGainAbs + 1.0f) + (fGainAbs - 1.0f) * fCos - fBeta * fSin);
            fA0 = (fGainAbs + 1.0f) - (fGainAbs - 1.0f) * fCos + fBeta * fSin;
            *out_pA1 = 2.0f * ((fGainAbs - 1.0f) - (fGainAbs + 1.0f) * fCos);
            *out_pA2 = (fGainAbs + 1.0f) - (fGainAbs - 1.0f) * fCos - fBeta * fSin;
            break;
        case AllPass:
            *out_pB0 = 1.0f - fAlpha;
            *out_pB1 = -2.0f * fCos;
            *out_pB2 = 1.0f + fAlpha;
            fA0 = 1.0f + fAlpha;
            *out_pA1 = -2.0f * fCos;
            *out_pA2 = 1.0f - fAlpha;
            break;
        }
        *out_pB0 /= fA0;
        *out_pB1 /= fA0;
        *out_pB2 /= fA0;
        *out_pA1 /= fA0;
        *out_pA2 /= fA0;
    }

    /// Based on heavy static: https://github.com/Wasted-Audio/hvcc/tree/develop/hvcc/generators/ir2c/static
    /// But only supports SIMD operations
    /// TODO: AVX support
    class Biquad
    {
    public:
        Biquad()
        {
        }

        void UpdateCoefficients(AkUInt32 in_uSampleRate,
                                FilterTypes in_eFilterType,
                                AkReal32 in_fFreq,
                                AkReal32 in_fQ,
                                AkReal32 in_fGainDB)
        {
            AkReal32 b0, b1, b2, a1, a2;
            ComputeBiquadCoefficients(in_uSampleRate,
                                      in_eFilterType,
                                      in_fFreq,
                                      in_fQ,
                                      in_fGainDB,
                                      &b0, &b1, &b2, &a1, &a2);
            UpdateCoefficients(b0, b1, b2, a1, a2);
        }

        void UpdateCoefficients(AkReal32 in_fB0, AkReal32 in_fB1, AkReal32 in_fB2, AkReal32 in_fA1, AkReal32 in_fA2)
        {
            const AkReal32 b0 = in_fB0;
            const AkReal32 b1 = in_fB1;
            const AkReal32 b2 = in_fB2;
            const AkReal32 a1 = -in_fA1;
            const AkReal32 a2 = -in_fA2;

            AkReal32 coeffs[4][8] =
            {
                {0, 0, 0, b0, b1, b2, a1, a2},
                {0, 0, b0, b1, b2, 0, a2, 0},
                {0, b0, b1, b2, 0, 0, 0, 0},
                {b0, b1, b2, 0, 0, 0, 0, 0},
            };

            for (int i = 0; i < 8; i++)
            {
                coeffs[1][i] += a1 * coeffs[0][i];
                coeffs[2][i] += a1 * coeffs[1][i] + a2 * coeffs[0][i];
                coeffs[3][i] += a1 * coeffs[2][i] + a2 * coeffs[1][i];
            }

            coeff_xp3 = AKSIMD_SETV_V4F32(coeffs[3][0], coeffs[2][0], coeffs[1][0], coeffs[0][0]);
            coeff_xp2 = AKSIMD_SETV_V4F32(coeffs[3][1], coeffs[2][1], coeffs[1][1], coeffs[0][1]);
            coeff_xp1 = AKSIMD_SETV_V4F32(coeffs[3][2], coeffs[2][2], coeffs[1][2], coeffs[0][2]);
            coeff_x0 = AKSIMD_SETV_V4F32(coeffs[3][3], coeffs[2][3], coeffs[1][3], coeffs[0][3]);
            coeff_xm1 = AKSIMD_SETV_V4F32(coeffs[3][4], coeffs[2][4], coeffs[1][4], coeffs[0][4]);
            coeff_xm2 = AKSIMD_SETV_V4F32(coeffs[3][5], coeffs[2][5], coeffs[1][5], coeffs[0][5]);
            coeff_ym1 = AKSIMD_SETV_V4F32(coeffs[3][6], coeffs[2][6], coeffs[1][6], coeffs[0][6]);
            coeff_ym2 = AKSIMD_SETV_V4F32(coeffs[3][7], coeffs[2][7], coeffs[1][7], coeffs[0][7]);
        }

        void CopyCoefficients(const Biquad* in_pOther)
        {
            coeff_xp3 = in_pOther->coeff_xp3;
            coeff_xp2 = in_pOther->coeff_xp2;
            coeff_xp1 = in_pOther->coeff_xp1;
            coeff_x0 = in_pOther->coeff_x0;
            coeff_xm1 = in_pOther->coeff_xm1;
            coeff_xm2 = in_pOther->coeff_xm2;
            coeff_ym1 = in_pOther->coeff_ym1;
            coeff_ym2 = in_pOther->coeff_ym2;
        }

        void Reset()
        {
            xm1 = AKSIMD_SETZERO_V4F32();
            xm2 = AKSIMD_SETZERO_V4F32();
            ym1 = AKSIMD_SETZERO_V4F32();
            ym2 = AKSIMD_SETZERO_V4F32();
        }

        AKSIMD_V4F32 Process(AKSIMD_V4F32 in_vIn)
        {
            AKSIMD_V4F32 x3 = AKSIMD_SHUFFLE_V4F32(in_vIn, in_vIn, AKSIMD_SHUFFLE(3,3,3,3));
            AKSIMD_V4F32 x2 = AKSIMD_SHUFFLE_V4F32(in_vIn, in_vIn, AKSIMD_SHUFFLE(2,2,2,2));
            AKSIMD_V4F32 x1 = AKSIMD_SHUFFLE_V4F32(in_vIn, in_vIn, AKSIMD_SHUFFLE(1,1,1,1));
            AKSIMD_V4F32 x0 = AKSIMD_SHUFFLE_V4F32(in_vIn, in_vIn, AKSIMD_SHUFFLE(0,0,0,0));

            AKSIMD_V4F32 a = AKSIMD_MUL_V4F32(coeff_xp3, x3);
            AKSIMD_V4F32 b = AKSIMD_MUL_V4F32(coeff_xp2, x2);
            AKSIMD_V4F32 c = AKSIMD_MUL_V4F32(coeff_xp1, x1);
            AKSIMD_V4F32 d = AKSIMD_MUL_V4F32(coeff_x0, x0);
            AKSIMD_V4F32 e = AKSIMD_MUL_V4F32(coeff_xm1, xm1);
            AKSIMD_V4F32 f = AKSIMD_MUL_V4F32(coeff_xm2, xm2);
            AKSIMD_V4F32 g = AKSIMD_MUL_V4F32(coeff_ym1, ym1);
            AKSIMD_V4F32 h = AKSIMD_MUL_V4F32(coeff_ym2, ym2);
            AKSIMD_V4F32 i = AKSIMD_ADD_V4F32(a, b);
            AKSIMD_V4F32 j = AKSIMD_ADD_V4F32(c, d);
            AKSIMD_V4F32 k = AKSIMD_ADD_V4F32(e, f);
            AKSIMD_V4F32 l = AKSIMD_ADD_V4F32(g, h);
            AKSIMD_V4F32 m = AKSIMD_ADD_V4F32(i, j);
            AKSIMD_V4F32 n = AKSIMD_ADD_V4F32(k, l);

            AKSIMD_V4F32 y = AKSIMD_ADD_V4F32(m, n);

            xm1 = x3;
            xm2 = x2;
            ym1 = AKSIMD_SHUFFLE_V4F32(y, y, AKSIMD_SHUFFLE(3,3,3,3));
            ym2 = AKSIMD_SHUFFLE_V4F32(y, y, AKSIMD_SHUFFLE(2,2,2,2));
            return y;
        }

        void Process(AkSampleType* io_pVector)
        {
            AKSIMD_V4F32 vIn = AKSIMD_LOAD_V4F32(io_pVector);
            AKSIMD_V4F32 vOut = Process(vIn);
            AKSIMD_STORE_V4F32(io_pVector, vOut);
        }

        void Process(AkSampleType* io_pBuffer, AkUInt32 in_uSize)
        {
            AKASSERT(in_uSize % uSimdSize == 0);

            for (AkUInt32 i = 0; i < in_uSize; i += uSimdSize)
            {
                Process(io_pBuffer + i);
            }
        }

        void Process(AkSampleType* in_pInBuffer, AkSampleType* out_pOutBuffer, AkUInt32 in_uSize)
        {
            AKASSERT(in_uSize % uSimdSize == 0);

            for (AkUInt32 i = 0; i < in_uSize; i += uSimdSize)
            {
                AKSIMD_V4F32 vIn = AKSIMD_LOAD_V4F32(in_pInBuffer + i);
                AKSIMD_V4F32 vOut = Process(vIn);
                AKSIMD_STORE_V4F32(out_pOutBuffer + i, vOut);
            }
        }

    private:
        static constexpr auto uSimdSize = 4;

        AKSIMD_V4F32 coeff_xp3{};
        AKSIMD_V4F32 coeff_xp2{};
        AKSIMD_V4F32 coeff_xp1{};
        AKSIMD_V4F32 coeff_x0{};
        AKSIMD_V4F32 coeff_xm1{};
        AKSIMD_V4F32 coeff_xm2{};
        AKSIMD_V4F32 coeff_ym1{};
        AKSIMD_V4F32 coeff_ym2{};

        // filter state
        AKSIMD_V4F32 xm1{};
        AKSIMD_V4F32 xm2{};
        AKSIMD_V4F32 ym1{};
        AKSIMD_V4F32 ym2{};
    };
}
