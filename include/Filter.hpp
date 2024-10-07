#pragma once
#include <vector>

#include "Biquad.hpp"


namespace Wpe
{
    class Filter
    {
    public:
        void Init(const AkUInt32 in_uSampleRate, const AkUInt32 in_uNumChannels)
        {
            m_uSampleRate = in_uSampleRate;
            m_biquads.resize(in_uNumChannels);
        }

        void Reset()
        {
            for (auto& biquad : m_biquads)
            {
                biquad.Reset();
            }
        }

        void Term()
        {
            m_biquads.clear();
        }

        void Resize(const AkUInt32 in_uNumChannels)
        {
            auto uOrgLength = m_biquads.size();
            m_biquads.resize(in_uNumChannels);
            if (in_uNumChannels > uOrgLength)
            {
                if (uOrgLength == 0)
                {
                    m_biquads[0].UpdateCoefficients(m_uSampleRate,
                                                    m_lastParams.eType,
                                                    m_lastParams.uFreq,
                                                    m_lastParams.fQ,
                                                    m_lastParams.fGainDB);
                    uOrgLength = 1;
                }
                for (AkUInt32 i = uOrgLength; i < in_uNumChannels; ++i)
                {
                    m_biquads[i].CopyCoefficients(&m_biquads[0]);
                }
            }
        }

        void Config(const FilterTypes in_eFilterType,
                    const AkReal32 in_fFreq,
                    const AkReal32 in_fQ,
                    const AkReal32 in_fGainDB)
        {
            m_lastParams.eType = in_eFilterType;
            m_lastParams.uFreq = in_fFreq;
            m_lastParams.fQ = in_fQ;
            m_lastParams.fGainDB = in_fGainDB;
            if (m_biquads.empty())
            {
                return;
            }
            auto it = m_biquads.begin();
            auto pFirstBiquad = &(*it);
            pFirstBiquad->UpdateCoefficients(m_uSampleRate, in_eFilterType, in_fFreq, in_fQ, in_fGainDB);
            ++it;
            for (; it != m_biquads.end(); ++it)
            {
                it->CopyCoefficients(pFirstBiquad);
            }
        }

        void ProcessChannel(AkSampleType* io_pBuffer, const AkUInt32 in_uChannelIndex, const AkUInt32 in_uSize)
        {
            if (AK_EXPECT_FALSE(in_uChannelIndex >= m_biquads.size()))
            {
                AKASSERT(false);
                return;
            }
            auto& biquad = m_biquads[in_uChannelIndex];
            biquad.Process(io_pBuffer, in_uSize);
        }

        void Process(AkAudioBuffer* io_pBuffer)
        {
            AKASSERT(io_pBuffer->NumChannels() == m_biquads.size());

            for (AkUInt32 i = 0; i < io_pBuffer->NumChannels(); ++i)
            {
                auto& biquad = m_biquads[i];
                biquad.Process(io_pBuffer->GetChannel(i), io_pBuffer->uValidFrames);
            }
        }

        void Process(AkAudioBuffer* in_pBuffer, AkAudioBuffer* out_pBuffer)
        {
            AKASSERT(out_pBuffer->NumChannels() == m_biquads.size());

            for (AkUInt32 i = 0; i < out_pBuffer->NumChannels(); ++i)
            {
                auto& biquad = m_biquads[i];
                biquad.Process(in_pBuffer->GetChannel(i), out_pBuffer->GetChannel(i), in_pBuffer->uValidFrames);
            }
        }

    private:
        struct FilterParams
        {
            FilterTypes eType;
            AkReal32 uFreq;
            AkReal32 fQ;
            AkReal32 fGainDB;
        };

        AkUInt32 m_uSampleRate = 0;
        FilterParams m_lastParams{};
        std::vector<Biquad> m_biquads;
    };
}
