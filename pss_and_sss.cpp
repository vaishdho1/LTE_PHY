/*
//  PSS_detection.cpp
//  LTE_PHY
The file contains the PSS and SSS generation and performs the following:
- Mapping to the symbols in the frequency domain.
- Converting into the time domain samples.
- PSS detection and symbol fine timing.
- SSS generationa and mapping into the frequency domain
- SSS detection and finding the half frame start index.
//Created by Vaishnavi  on 05/01/20.
//Copyright © 2020 none. All rights reserved.
*/


#include "liblte_phy.h"
#include<iostream>
#include "liblte_phy.h"
#include<math.h>
using namespace std;
void symbols_to_samples(LIBLTE_PHY *phy_struct,float *symb_re,float *symb_im,int symbol_offset,float *samps_re,float *samps_im, int *N_samps)
{
    int CP_len;
    int i,FFT_SIZE,FFT_padsize;
    FFT_SIZE=128;
    FFT_padsize= (FFT_SIZE-72)/2;

fftw_complex *s2s_in;
fftw_complex *s2s_out;

fftw_plan symbs_to_samps_dl_plan;
//fftw_plan samps_to_symbs_dl_plan;
s2s_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * phy_struct->N_samps_per_symb);
s2s_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * phy_struct->N_samps_per_symb);
symbs_to_samps_dl_plan= fftw_plan_dft_1d(phy_struct->N_samps_per_symb,s2s_in,s2s_out,FFTW_BACKWARD,FFTW_MEASURE);

if((symbol_offset%7)==0)
    CP_len=phy_struct->N_samps_cp_1_0;
else
    CP_len=phy_struct->N_samps_cp_1_else;

    for(i=0;i<phy_struct->N_samps_per_symb;i++)
    {
        s2s_in[i][0]=0;
        s2s_in[i][1]=0;
    }

    for(i=0;i<FFT_SIZE/2-FFT_padsize;i++)
    {   //Positive spectra
        s2s_in[i][0]= symb_re[i+FFT_SIZE/2-FFT_padsize];
        s2s_in[i][1]= symb_im[i+FFT_SIZE/2-FFT_padsize];

        //Negative spectra
        s2s_in[phy_struct->N_samps_per_symb-i-1][0] =symb_re[FFT_SIZE/2-FFT_padsize-i-1];
        s2s_in[phy_struct->N_samps_per_symb-i-1][1] =symb_im[FFT_SIZE/2-FFT_padsize-i-1];
    }


    fftw_execute(symbs_to_samps_dl_plan);

    for(i=0;i<phy_struct->N_samps_per_symb;i++)
    {
        samps_re[i+CP_len]=s2s_out[i][0];
        samps_im[i+CP_len]=s2s_out[i][1];

    }

    for(i=0;i<CP_len;i++)
    {
        samps_re[i]=samps_re[phy_struct->N_samps_per_symb+i];
        samps_im[i]=samps_im[phy_struct->N_samps_per_symb+i];
    }

    *N_samps =phy_struct->N_samps_per_symb+CP_len;
    cout<<"The total no. of samples of the symbol are "<<*N_samps<<endl;


}
void samples_to_symbols(LIBLTE_PHY *phy_struct,float *symb_re,float *symb_im,int symbol_offset,int slot_start_idx,float *samps_re,float *samps_im)
{
    int CP_len,index;
    int i,FFT_SIZE,FFT_padsize;
    FFT_SIZE=128;
    FFT_padsize= (FFT_SIZE-72)/2;

fftw_complex *s2s_in;
fftw_complex *s2s_out;

fftw_plan samps_to_symbs_dl_plan;
//fftw_plan samps_to_symbs_dl_plan;
s2s_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * phy_struct->N_samps_per_symb);
s2s_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * phy_struct->N_samps_per_symb);
samps_to_symbs_dl_plan= fftw_plan_dft_1d(phy_struct->N_samps_per_symb,s2s_in,s2s_out,FFTW_FORWARD,FFTW_MEASURE);

if((symbol_offset%7)==0)
    CP_len=phy_struct->N_samps_cp_1_0;
else
    CP_len=phy_struct->N_samps_cp_1_else;
    cout<<"CP_len"<<CP_len<<endl;
    index=slot_start_idx+(128+9)*symbol_offset;
    for(i=0;i<phy_struct->N_samps_per_symb;i++)
    {
        s2s_in[i][0]=samps_re[index+i+CP_len-1];
        s2s_in[i][1]=samps_im[index+i+CP_len-1];

    }

    fftw_execute(samps_to_symbs_dl_plan);


    for(i=0;i<FFT_SIZE/2-FFT_padsize;i++)
    {   //Positive spectra
        symb_re[i+FFT_SIZE/2-FFT_padsize]=s2s_out[i][0];
        symb_im[i+FFT_SIZE/2-FFT_padsize]=s2s_out[i][1];

        //Negative spectra
         symb_re[FFT_SIZE/2-FFT_padsize-i-1]=s2s_out[phy_struct->N_samps_per_symb-i-1][0];
        symb_im[FFT_SIZE/2-FFT_padsize-i-1]= s2s_out[phy_struct->N_samps_per_symb-i-1][1];
    }


}
void pss_generate(int N_id2, float* pss_re,float* pss_im)
{
    float root;
    int i;
    if(N_id2==0)
      root=25;
else if(N_id2==1)
    root=29;
else root=34;
    for(i=0;i<31;i++)
    { pss_re[i]=cosf(-M_PI*root*i*(i+1)/63);
      pss_im[i]=sinf(-M_PI*root*i*(i+1)/63);
    }
    for(i=31;i<62;i++)
    { pss_re[i]=cosf(-M_PI*root*(i+2)*(i+1)/63);
      pss_im[i]=sinf(-M_PI*root*(i+2)*(i+1)/63);
    }

    cout<<"The pss generation is successful"<<endl;
}


void find_pss_and_fine_timing(LIBLTE_PHY *phy_struct,
                                                      float             *i_samps,
                                                      float             *q_samps,
                                                      int            *symb_starts,
                                                      int            *N_id_2,
                                                      int            *pss_symb,
                                                      float             *pss_thresh,
                                                      float             *freq_offset)
{

    float              corr_re_n1;
    float              corr_im_n1;
    float              corr_re;
    float              corr_im;
    float              corr_re_p1;
    float              corr_im_p1;
    float              abs_corr_n1;
    float              abs_corr;
    float              abs_corr_p1;
    float              corr_max;
    float              pss_re[63];
    float              pss_im[63];
    float             *pss_mod_re;
    float             *pss_mod_im;
    int              i;
    int              idx = 0;
    int             j;
    int             k;
    int             z;
    int             N_s;
    int             N_symb;
    int             pss_timing_idx;
    int               timing;
    int            N_samples_per_symbol;

      // Generate PSS
        for(i=0; i<3; i++)
        {
            for(j=0; j<72; j++)
            {
                phy_struct->pss_mod_re[i][j] = 0;
                phy_struct->pss_mod_im[i][j] = 0;
            }
            pss_generate(i, pss_re, pss_im);
            for(j=0; j<62; j++)
            {
                k                                 = j - 31 + 36;
                phy_struct->pss_mod_re_n1[i][k-1] = pss_re[j];
                phy_struct->pss_mod_im_n1[i][k-1] = pss_im[j];
                phy_struct->pss_mod_re[i][k]      = pss_re[j];
                phy_struct->pss_mod_im[i][k]      = pss_im[j];
                phy_struct->pss_mod_re_p1[i][k+1] = pss_re[j];
                phy_struct->pss_mod_im_p1[i][k+1] = pss_im[j];
            }
        }

        // Demod symbols and correlate with PSS
    corr_max = 0;int l;
        for(i=0; i<12; i++)
        {
            for(j=0; j<7; j++)
            { if(j%7==0)
                N_samples_per_symbol=138;
              else
                  N_samples_per_symbol=137;

                samples_to_symbols (phy_struct,
                                    phy_struct->rx_symb_re,
                                    phy_struct->rx_symb_im,
                                    0,
                                    symb_starts[j]+(phy_struct->N_samps_per_slot*i),
                                    &i_samps[0],
                                    &q_samps[0]
                                    );

                for(k=0; k<3; k++)
                {
                    corr_re_n1 = 0;
                    corr_im_n1 = 0;
                    corr_re    = 0;
                    corr_im    = 0;
                    corr_re_p1 = 0;
                    corr_im_p1 = 0;
                    for(z=0; z<72; z++)
                    {
                        corr_re_n1 += (phy_struct->rx_symb_re[z]*phy_struct->pss_mod_re_n1[k][z] +
                                       phy_struct->rx_symb_im[z]*phy_struct->pss_mod_im_n1[k][z]);
                        corr_im_n1 += (phy_struct->rx_symb_re[z]*phy_struct->pss_mod_im_n1[k][z] -
                                       phy_struct->rx_symb_im[z]*phy_struct->pss_mod_re_n1[k][z]);
                        corr_re    += (phy_struct->rx_symb_re[z]*phy_struct->pss_mod_re[k][z] +
                                       phy_struct->rx_symb_im[z]*phy_struct->pss_mod_im[k][z]);
                        corr_im    += (phy_struct->rx_symb_re[z]*phy_struct->pss_mod_im[k][z] -
                                       phy_struct->rx_symb_im[z]*phy_struct->pss_mod_re[k][z]);
                        corr_re_p1 += (phy_struct->rx_symb_re[z]*phy_struct->pss_mod_re_p1[k][z] +
                                       phy_struct->rx_symb_im[z]*phy_struct->pss_mod_im_p1[k][z]);
                        corr_im_p1 += (phy_struct->rx_symb_re[z]*phy_struct->pss_mod_im_p1[k][z] -
                                       phy_struct->rx_symb_im[z]*phy_struct->pss_mod_re_p1[k][z]);
                    }
                    abs_corr_n1 = sqrt(corr_re_n1*corr_re_n1 + corr_im_n1*corr_im_n1);
                    abs_corr    = sqrt(corr_re*corr_re + corr_im*corr_im);
                    abs_corr_p1 = sqrt(corr_re_p1*corr_re_p1 + corr_im_p1*corr_im_p1);
                    if(abs_corr_n1 > corr_max)
                    {
                        idx       = -1;
                        corr_max  = abs_corr_n1;
                        *pss_symb = (i*7)+j;
                        *N_id_2   = k;
                    }
                    if(abs_corr > corr_max)
                    {
                        idx       = 0;
                        corr_max  = abs_corr;
                        *pss_symb = (i*7)+j;
                        *N_id_2   = k;
                    }
                    if(abs_corr_p1 > corr_max)
                    {
                        idx       = 1;
                        corr_max  = abs_corr_p1;
                        *pss_symb = (i*7)+j;
                        *N_id_2   = k;
                    }
                }
            }
        }
        if(-1 == idx)
        {
            pss_mod_re   = &phy_struct->pss_mod_re_n1[*N_id_2][0];
            pss_mod_im   = &phy_struct->pss_mod_im_n1[*N_id_2][0];
            *freq_offset = -15000; // FIXME
        }else if(0 == idx){
            pss_mod_re   = &phy_struct->pss_mod_re[*N_id_2][0];
            pss_mod_im   = &phy_struct->pss_mod_im[*N_id_2][0];
            *freq_offset = 0;
        }else{
            pss_mod_re   = &phy_struct->pss_mod_re_p1[*N_id_2][0];
            pss_mod_im   = &phy_struct->pss_mod_im_p1[*N_id_2][0];
            *freq_offset = 15000; // FIXME
        }


        // Find optimal timing
        N_s      = (*pss_symb)/7;
        N_symb   = (*pss_symb)%7;
        corr_max = 0;
        timing   = 0;
        for(i=-40; i<40; i++)
        {
            idx = symb_starts[N_symb] + (phy_struct->N_samps_per_slot*N_s);
            if(i < 0)
            {
                if(idx >= -i)
                {
                    idx += i;
                }
            }else{
                idx += i;
            }

                         samples_to_symbols (phy_struct,
                                              phy_struct->rx_symb_re,
                                              phy_struct->rx_symb_im,
                                              0,
                                              idx,
                                              &i_samps[0],
                                              &q_samps[0]
                                              );
            corr_re = 0;
            corr_im = 0;
            for(j=0; j<72; j++)
            {
                corr_re += (phy_struct->rx_symb_re[j]*pss_mod_re[j] +
                            phy_struct->rx_symb_im[j]*pss_mod_im[j]);
                corr_im += (phy_struct->rx_symb_re[j]*pss_mod_im[j] -
                            phy_struct->rx_symb_im[j]*pss_mod_re[j]);
            }
            abs_corr = sqrt(corr_re*corr_re + corr_im*corr_im);
            if(abs_corr > corr_max)
            {
                corr_max = abs_corr;
                timing   = i;
            }
        }

        *pss_thresh = corr_max;
    cout<<"pss_thresh"<<*pss_thresh<<endl;

        // Construct fine symbol start locations
        pss_timing_idx = symb_starts[N_symb]+(phy_struct->N_samps_per_slot*N_s)+timing;

        while((pss_timing_idx + phy_struct->N_samps_per_symb + phy_struct->N_samps_cp_l_else) < phy_struct->N_samps_per_slot)
        {
            pss_timing_idx += phy_struct->N_samps_per_frame;
        }

        symb_starts[0] = pss_timing_idx + (phy_struct->N_samps_per_symb+phy_struct->N_samps_cp_l_else)*1 - phy_struct->N_samps_per_slot;
    for(i=1;i<7;i++)
        {
            symb_starts[i] = pss_timing_idx + (phy_struct->N_samps_per_symb+phy_struct->N_samps_cp_l_else)*i + phy_struct->N_samps_per_symb+phy_struct->N_samps_cp_l_0 - phy_struct->N_samps_per_slot;
        }
    cout<<"PSS timing"<<endl;
    for(i=0;i<7;i++)
        cout<<"symb_starts"<<symb_starts[i]<<endl;
    cout<<"N_id_2"<<*N_id_2<<endl;
    }

void generate_sss(LIBLTE_PHY *phy_struct,int N_id_1,int N_id_2,float *sss_re_0,float *sss_im_0,float *sss_re_5,float *sss_im_5)
{
    int q_prime,m_prime,m_0,m_1,i,q;

    q_prime=N_id_1/30;
    q=(N_id_1+(q_prime*(q_prime+1)/2))/30;
    m_prime=N_id_1+q*(q+1)/2;
    m_0=m_prime%31;
    m_1=(m_0+(m_prime/31)+1)%31;
    memset(phy_struct->sss_x_s_tilda,0,31);
    for(i=0;i<31;i++)

    phy_struct->sss_x_s_tilda[4]=1;
    for(i=0;i<26;i++)
        phy_struct->sss_x_s_tilda[i+5]=(phy_struct->sss_x_s_tilda[i+2]+phy_struct->sss_x_s_tilda[i])%2;

    for(i=0;i<31;i++)
        phy_struct->sss_s_tilda[i]=1-2*(phy_struct->sss_x_s_tilda[i]);



    memset(phy_struct->sss_x_c_tilda,0,31);
    phy_struct->sss_x_c_tilda[4]=1;

    for(i=0;i<26;i++)
        phy_struct->sss_x_c_tilda[i+5]=(phy_struct->sss_x_c_tilda[i+3]+phy_struct->sss_x_c_tilda[i])%2;
    for(i=0;i<31;i++)
        phy_struct->sss_c_tilda[i]=1-2*(phy_struct->sss_x_c_tilda[i]);



    memset(phy_struct->sss_x_z_tilda,0,31);
       phy_struct->sss_x_z_tilda[4]=1;

       for(i=0;i<26;i++)
       phy_struct->sss_x_z_tilda[i+5]=(phy_struct->sss_x_z_tilda[i+4]+phy_struct->sss_x_z_tilda[i+2]+phy_struct->sss_x_z_tilda[i+1]+phy_struct->sss_x_z_tilda[i])%2;

       for(i=0;i<31;i++)
           phy_struct->sss_z_tilda[i]=1-2*(phy_struct->sss_x_z_tilda[i]);

    for(i=0;i<31;i++)
    {
        phy_struct->sss_s0_m0[i]=phy_struct->sss_s_tilda[(i+m_0)%31];
        phy_struct->sss_s1_m1[i]=phy_struct->sss_s_tilda[(i+m_1)%31];
    }

    for(i=0;i<31;i++)
    {
        phy_struct->sss_c0[i]=phy_struct->sss_c_tilda[(i+N_id_2)%31];
        phy_struct->sss_c1[i]=phy_struct->sss_c_tilda[(i+N_id_2+3)%31];
    }

    for(i=0;i<31;i++)
    {
        phy_struct->sss_z1_m0[i]=phy_struct->sss_z_tilda[(i+(m_0%8))%31];
        phy_struct->sss_z1_m1[i]=phy_struct->sss_z_tilda[(i+(m_1%8))%31];
    }

    for(i=0;i<31;i++)
    {
        sss_re_0[2*i]=phy_struct->sss_s0_m0[i]*phy_struct->sss_c0[i];
        sss_re_5[2*i]=phy_struct->sss_s1_m1[i]*phy_struct->sss_c0[i];
        sss_im_0[2*i]=0;
        sss_im_5[2*i]=0;
        sss_im_0[2*i+1]=0;
        sss_im_5[2*i+1]=0;
        sss_re_0[2*i+1]=phy_struct->sss_s1_m1[i]*phy_struct->sss_c1[i]*phy_struct->sss_z1_m0[i];
        sss_re_5[2*i+1]=phy_struct->sss_s0_m0[i]*phy_struct->sss_c1[i]*phy_struct->sss_z1_m1[i];

    }





}


void sss_mapping(LIBLTE_PHY *phy_struct, FRAME *subframe, int N_id_1, int N_id_2,float *samps_re,float *samps_im)
{ int k,i,N_samps;
    generate_sss(phy_struct,N_id_1,N_id_2,&phy_struct->sss_mod_re_0[0][0],&phy_struct->sss_mod_im_0[0][0],&phy_struct->sss_mod_re_5[0][0],&phy_struct->sss_mod_im_5[0][0]);
    subframe->num=5;

    if(subframe->num==0)
    {cout<<"Entered 1"<<endl;
        for(i=0;i<62;i++)
        {
            k=i-31+36;

            subframe->tx_symb_re[5][k]=phy_struct->sss_mod_re_0[0][i];
            subframe->tx_symb_im[5][k]=phy_struct->sss_mod_im_0[0][i];

        }

        }
else if (subframe->num==5)
{cout<<"Entered 5"<<endl;
    for(i=0;i<62;i++)
    {
        k=i-31+36;

        subframe->tx_symb_re[75][k]=phy_struct->sss_mod_re_5[0][i];
        subframe->tx_symb_im[5][k]=phy_struct->sss_mod_im_5[0][i];

    }

}
   symbols_to_samples(phy_struct,&subframe->tx_symb_re[75][0],&subframe->tx_symb_im[75][0],5,&samps_re[0],&samps_im[0],&N_samps);

    cout<<"sss mapping successful"<<endl;
}

void find_sss(LIBLTE_PHY *phy_struct,
                                      float             *i_samps,
                                      float             *q_samps,
                                      int             N_id_2,
                                      int            *symb_starts,
                                      float              pss_thresh,
                                      int            *N_id_1,
                                      int            *frame_start_idx)
{

    float             sss_thresh;
    float             corr_re;
    float             corr_im;
    float             abs_corr;
    int            i;
    int            j;
    int            k;


        // Generate secondary synchronization signals
        memset(phy_struct->sss_mod_re_0, 0, sizeof(float)*168*6*12);
        memset(phy_struct->sss_mod_im_0, 0, sizeof(float)*168*6*12);
        memset(phy_struct->sss_mod_re_5, 0, sizeof(float)*168*6*12);
        memset(phy_struct->sss_mod_im_5, 0, sizeof(float)*168*6*12);
        for(i=0; i<168; i++)
        {
            generate_sss(phy_struct,
                         i,
                         N_id_2,
                         phy_struct->sss_re_0,
                         phy_struct->sss_im_0,
                         phy_struct->sss_re_5,
                         phy_struct->sss_im_5);
            for(j=0; j<62; j++)
            {
                k                              = j - 31 + (72)/2;
                phy_struct->sss_mod_re_0[i][k] = phy_struct->sss_re_0[j];
                phy_struct->sss_mod_im_0[i][k] = phy_struct->sss_im_0[j];
                phy_struct->sss_mod_re_5[i][k] = phy_struct->sss_re_5[j];
                phy_struct->sss_mod_im_5[i][k] = phy_struct->sss_im_5[j];
            }
        }
        sss_thresh = pss_thresh * 0.9;
    cout<<"symb_starts[5]"<<symb_starts[5]<<endl;
        // Demod symbol and search for secondary synchronization signals
        samples_to_symbols(phy_struct,
                           phy_struct->rx_symb_re,
                           phy_struct->rx_symb_im,
                              0,
                              symb_starts[5],
                              &i_samps[0],
                              &q_samps[0]
                              );

        for(i=0; i<168; i++)
        {
            corr_re = 0;
            corr_im = 0;
            for(j=0; j<72; j++)
            {
                corr_re += (phy_struct->rx_symb_re[j]*phy_struct->sss_mod_re_0[i][j] +
                            phy_struct->rx_symb_im[j]*phy_struct->sss_mod_im_0[i][j]);
                corr_im += (phy_struct->rx_symb_re[j]*phy_struct->sss_mod_im_0[i][j] -
                            phy_struct->rx_symb_im[j]*phy_struct->sss_mod_re_0[i][j]);
            }
            abs_corr = sqrt(corr_re*corr_re + corr_im*corr_im);

            if(abs_corr > sss_thresh)
            {
                while(symb_starts[5] < ((phy_struct->N_samps_per_symb + phy_struct->N_samps_cp_l_else)*4 + phy_struct->N_samps_per_symb + phy_struct->N_samps_cp_l_0))
                {
                    symb_starts[5] += phy_struct->N_samps_per_frame;
                }
                *N_id_1          = i;
                *frame_start_idx = symb_starts[5] - ((phy_struct->N_samps_per_symb + phy_struct->N_samps_cp_l_else)*4 + phy_struct->N_samps_per_symb + phy_struct->N_samps_cp_l_0);

                break;
            }

            corr_re = 0;
            corr_im = 0;
            for(j=0; j<72; j++)
            {
                corr_re += (phy_struct->rx_symb_re[j]*phy_struct->sss_mod_re_5[i][j] +
                            phy_struct->rx_symb_im[j]*phy_struct->sss_mod_im_5[i][j]);
                corr_im += (phy_struct->rx_symb_re[j]*phy_struct->sss_mod_im_5[i][j] -
                            phy_struct->rx_symb_im[j]*phy_struct->sss_mod_re_5[i][j]);
            }
            abs_corr = sqrt(corr_re*corr_re + corr_im*corr_im);

            if(abs_corr > sss_thresh)
            {
                while(symb_starts[5] < (((phy_struct->N_samps_per_symb + phy_struct->N_samps_cp_l_else)*4 + phy_struct->N_samps_per_symb + phy_struct->N_samps_cp_l_0) + phy_struct->N_samps_per_slot*10))
                {
                    symb_starts[5] += phy_struct->N_samps_per_frame;
                }
                *N_id_1          = i;
                *frame_start_idx = symb_starts[5] - ((phy_struct->N_samps_per_symb + phy_struct->N_samps_cp_l_else)*4 + phy_struct->N_samps_per_symb + phy_struct->N_samps_cp_l_0) - phy_struct->N_samps_per_slot*10;
                break;
            }

        }
    cout<<"The N_id_1 and frame start"<<endl;
    cout<<*N_id_1<<endl;
    cout<<*frame_start_idx<<endl;
    cout<<"sss finding successful"<<endl;

    }

int main()

{

    LIBLTE_PHY *phy_struct = new LIBLTE_PHY;
    float pss_re[63],pss_im[63];
    int pss_symb;
    float pss_thresh,freq_offset;
    int frame_start_index;
    int N_ID2=0;
    int N_ID1=34;
    int N_id1,N_id2;
    int i,k;
    int symb_starts[7]={0,138,275,412,549,686,823};
    float samps_re[137]={0};
    float isamps[11520]={0};
    float qsamps[11520]={0};
    float samps_im[137]={0};
    //float samps_re[137],samps_im[137];
    int *N_samps = new int;

         FRAME *subframe= new FRAME;

         pss_generate(N_ID2,pss_re,pss_im);
         //symbols_to_samples

         for(i=0;i<62;i++)
         {   k=i-31+36;
             subframe->tx_symb_re[76][k]=pss_re[i];
             subframe->tx_symb_im[76][k]=pss_im[i];
         }
    cout<<"Successful mapping of PSS symbols in frequency domain"<<endl;

    symbols_to_samples(phy_struct,&subframe->tx_symb_re[76][0],&subframe->tx_symb_im[76][0],6,samps_re,samps_im,N_samps);

    for(i=10423;i<10560;i++)
       {
           isamps[i]=samps_re[i-10423];
           qsamps[i]=samps_im[i-10423];
       }
    find_pss_and_fine_timing(phy_struct,&isamps[0],&qsamps[0],symb_starts,&N_id2,&pss_symb,&pss_thresh,&freq_offset);

    cout<<"The pss and fine timing detected"<<endl;

    cout<<"The samples in time domain for pss are successfully converted"<<endl;

    sss_mapping(phy_struct,subframe,N_ID1, N_ID2,&samps_re[0],&samps_im[0]);

    for(i=10286;i<10423;i++)
    {
        isamps[i]=samps_re[i-10286];
        qsamps[i]=samps_im[i-10286];

    }

    cout<<"The sss signals are mapped"<<endl;


        find_sss(phy_struct,&isamps[0],
             &qsamps[0],N_id2,&symb_starts[0],pss_thresh,&N_id1,&frame_start_index);
    cout<<"sss signal detection successful"<<endl;

}


