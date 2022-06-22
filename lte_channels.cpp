/*----------------------------------------------------------

  lte_channels.cpp
  LTE_PHY
  The file contains the implementations of all the downlink channels
  and their respective functions.
  The channel related functions are defined in
  Created by Vaishnavi  on 28/02/20.
  Copyright © 2020 none. All rights reserved.
Adding all the downlink phy channels processing and their mapping.
Improvements and bug fixes:
18/03 - Changed the pbch_channel_map function to do encoding only when SFN%4=0


 -----------------------------------------------------------------*/
#include "includes.h"
using namespace std;
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
    index=slot_start_idx+(phy_struct->N_samps_per_symb+phy_struct->N_samps_cp_1_else)*symbol_offset;
    if(symbol_offset>0)
    {
        index+=phy_struct->N_samps_cp_1_0-phy_struct->N_samps_cp_1_else;
    }
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


}
void conv_encode(LIBLTE_PHY *phy_struct,
                 int             *c_bits,
                 int             N_c_bits,
                 int             constraint_len,
                 int             rate,
                 int            *g,
                 bool               tail_bit,
                 int             *d_bits,
                 int            *N_d_bits)
{
    int i;
    int j;
    int k;
    int  s_reg[constraint_len];
    int  g_array[3][constraint_len];

    // Initialize the shift register
    if(tail_bit)
    {
        for(i=0; i<constraint_len; i++)
        {
            s_reg[i] = c_bits[N_c_bits-i-1];
        }
    }else{
        for(i=0; i<constraint_len; i++)
        {
            s_reg[i] = 0;
        }
    }

    // Convert g from octal to binary array
    for(i=0; i<rate; i++)
    {
        for(j=0; j<constraint_len; j++)
        {
            g_array[i][j] = (g[i] >> (constraint_len-j-1)) & 1;
        }
    }

    // Convolutionally encode input
    for(i=0; i<N_c_bits; i++)
    {
        // Add next bit to shift register
        for(j=constraint_len-1; j>0; j--)
        {
            s_reg[j] = s_reg[j-1];
        }
        s_reg[0] = c_bits[i];

        // Determine the output bits
        for(j=0; j<rate; j++)
        {
            d_bits[i*rate + j] = 0;

            for(k=0; k<constraint_len; k++)
            {
                d_bits[i*rate + j] += s_reg[k]*g_array[j][k];
            }
            d_bits[i*rate + j] %= 2;
        }
    }

    *N_d_bits = N_c_bits*rate;
}
void convert_to_bits(int val,int **bits,int N_bits)
{
    int i;

    for(i=0;i<N_bits;i++)
    {
        (*bits)[i]=((N_bits-i)<<val)&0x1;
    }
    *bits+=N_bits;
}

void scramb_generate(int  c_init,
int len,
int *c)
{
    int i;
    int x1;
    int x2;

    int  new_bit1;
    int  new_bit2;

    // Initialize the 2nd m-sequence
    x2 = c_init;

    // Advance the 2d m-sequence
    for(i=0; i<(1600-31); i++)
    {
        new_bit2 = ((x2 >> 3) ^ (x2 >> 2) ^ (x2 >> 1) ^ x2) & 0x1;

        x2 = (x2 >> 1) | (new_bit2 << 30);
    }

    // Initialize the 1st m-sequence
    x1 = 0x54D21B24; // This is the result of advancing the initial value of 0x00000001

    // Generate c
    for(i=0; i<len; i++)
    {
        new_bit1 = ((x1 >> 3) ^ x1) & 0x1;
        new_bit2 = ((x2 >> 3) ^ (x2 >> 2) ^ (x2 >> 1) ^ x2) & 0x1;

        x1 = (x1 >> 1) | (new_bit1 << 30);
        x2 = (x2 >> 1) | (new_bit2 << 30);

        c[i] = new_bit1 ^ new_bit2;
    }



}
void cfi_encode(LIBLTE_PHY *phy_struct,int cfi,int *bits,int *N_bits)
{
    int *cfi_bits;
    int i;
    *N_bits=32;
    if(cfi==1)
        cfi_bits=CFI_0;
    else if(cfi==2)
        cfi_bits=CFI_1;
    else if(cfi==3)
        cfi_bits=CFI_2;
    else
        cfi_bits=CFI_3;

    for(i=0;i<*N_bits;i++)
        bits[i]=cfi_bits[i];

}

void pre_coder_dl(float                          *x_re,
                  float                          *x_im,
                  int                          M_layer_symb,
                  int                           N_ant,
                  LIBLTE_PHY_PRE_CODER_TYPE_ENUM  type,
                  float                          *y_re,
                  float                          *y_im,
                  int                         y_len,
                  int                       *M_ap_symb)
{
    float  *x_re_ptr[N_ant];
    float  *x_im_ptr[N_ant];
    float  *y_re_ptr[N_ant];
    float  *y_im_ptr[N_ant];
    float   one_over_sqrt_2 = 1/sqrt(2);
    int32_t  i;
    int32_t  p;

    // Index all arrays   Layer mapping
    for(p=0; p<N_ant; p++)
    {
        x_re_ptr[p] = &x_re[p*M_layer_symb];
        x_im_ptr[p] = &x_im[p*M_layer_symb];
        y_re_ptr[p] = &y_re[p*y_len];
        y_im_ptr[p] = &y_im[p*y_len];
    }

    if(N_ant == 1)
    {
        // 3GPP TS 36.211 v10.1.0 section 6.3.4.1
        *M_ap_symb = M_layer_symb;
        for(i=0; i<*M_ap_symb; i++)
        {
            y_re_ptr[0][i] = x_re_ptr[0][i];
            y_im_ptr[0][i] = x_im_ptr[0][i];
        }
    }else if(N_ant == 2){
        if(type==1)
        {
            // 3GPP TS 36.211 v10.1.0 section 6.3.4.3
            *M_ap_symb = 2*M_layer_symb;
            for(i=0; i<M_layer_symb; i++)
            {
                y_re_ptr[0][2*i+0] = +one_over_sqrt_2 * x_re_ptr[0][i];
                y_im_ptr[0][2*i+0] = +one_over_sqrt_2 * x_im_ptr[0][i];
                y_re_ptr[1][2*i+0] = -one_over_sqrt_2 * x_re_ptr[1][i];
                y_im_ptr[1][2*i+0] = +one_over_sqrt_2 * x_im_ptr[1][i];
                y_re_ptr[0][2*i+1] = +one_over_sqrt_2 * x_re_ptr[1][i];
                y_im_ptr[0][2*i+1] = +one_over_sqrt_2 * x_im_ptr[1][i];
                y_re_ptr[1][2*i+1] = +one_over_sqrt_2 * x_re_ptr[0][i];
                y_im_ptr[1][2*i+1] = -one_over_sqrt_2 * x_im_ptr[0][i];
            }
        }else{
            // : Currently only supporting TX Diversity
        }
    }else{
        if(type==1)
        {
            // 3GPP TS 36.211 v10.1.0 section 6.3.4.3
            if(x_re_ptr[2][M_layer_symb-1] == TX_NULL_SYMB &&
               x_im_ptr[2][M_layer_symb-1] == TX_NULL_SYMB &&
               x_re_ptr[3][M_layer_symb-1] == TX_NULL_SYMB &&
               x_im_ptr[3][M_layer_symb-1] == TX_NULL_SYMB)
            {
                *M_ap_symb = 4*M_layer_symb - 2;
            }else{
                *M_ap_symb = 4*M_layer_symb;
            }
            for(i=0; i<M_layer_symb; i++)
            {
                y_re_ptr[0][4*i+0] = +one_over_sqrt_2 * x_re_ptr[0][i];
                y_im_ptr[0][4*i+0] = +one_over_sqrt_2 * x_im_ptr[0][i];
                y_re_ptr[1][4*i+0] = 0;
                y_im_ptr[1][4*i+0] = 0;
                y_re_ptr[2][4*i+0] = -one_over_sqrt_2 * x_re_ptr[1][i];
                y_im_ptr[2][4*i+0] = +one_over_sqrt_2 * x_im_ptr[1][i];
                y_re_ptr[3][4*i+0] = 0;
                y_im_ptr[3][4*i+0] = 0;
                y_re_ptr[0][4*i+1] = +one_over_sqrt_2 * x_re_ptr[1][i];
                y_im_ptr[0][4*i+1] = +one_over_sqrt_2 * x_im_ptr[1][i];
                y_re_ptr[1][4*i+1] = 0;
                y_im_ptr[1][4*i+1] = 0;
                y_re_ptr[2][4*i+1] = +one_over_sqrt_2 * x_re_ptr[0][i];
                y_im_ptr[2][4*i+1] = -one_over_sqrt_2 * x_im_ptr[0][i];
                y_re_ptr[3][4*i+1] = 0;
                y_im_ptr[3][4*i+1] = 0;
                y_re_ptr[0][4*i+2] = 0;
                y_im_ptr[0][4*i+2] = 0;
                y_re_ptr[1][4*i+2] = +one_over_sqrt_2 * x_re_ptr[2][i];
                y_im_ptr[1][4*i+2] = +one_over_sqrt_2 * x_im_ptr[2][i];
                y_re_ptr[2][4*i+2] = 0;
                y_im_ptr[2][4*i+2] = 0;
                y_re_ptr[3][4*i+2] = -one_over_sqrt_2 * x_re_ptr[3][i];
                y_im_ptr[3][4*i+2] = +one_over_sqrt_2 * x_im_ptr[3][i];
                y_re_ptr[0][4*i+3] = 0;
                y_im_ptr[0][4*i+3] = 0;
                y_re_ptr[1][4*i+3] = +one_over_sqrt_2 * x_re_ptr[3][i];
                y_im_ptr[1][4*i+3] = +one_over_sqrt_2 * x_im_ptr[3][i];
                y_re_ptr[2][4*i+3] = 0;
                y_im_ptr[2][4*i+3] = 0;
                y_re_ptr[3][4*i+3] = +one_over_sqrt_2 * x_re_ptr[2][i];
                y_im_ptr[3][4*i+3] = -one_over_sqrt_2 * x_im_ptr[2][i];
            }
        }else{
            // FIXME: Currently only supporting TX Diversity
        }
    }
}

void layer_mapper(float        *d_re,
float                          *d_im,
int                          M_symb,
int                           N_ant,
int                          N_codewords,
LIBLTE_PHY_PRE_CODER_TYPE_ENUM  type,
float                          *x_re,
float                          *x_im,
int                         *M_layer_symb)
{
    int i;
    if(N_ant==1&&N_codewords==1)
    {
        *M_layer_symb=M_symb;
        for(i=0;i<*M_layer_symb;i++)
        {
            x_re[i]=d_re[i];
        x_im[i]=d_im[i];
        }
    }
    else if(type==1)
    {
        if(N_ant==2)
        {*M_layer_symb=M_symb/2 ;
            for(i=0;i<*M_layer_symb;i++)
            {
                x_re[i]=d_re[2*i];
                x_im[i]=d_im[2*i];
                x_re[i+*M_layer_symb]=d_re[2*i+1];
                x_im[i+*M_layer_symb]=d_im[2*i+1];
            }
        }
        else if(N_ant==4)
        {if(M_symb%4==0)
            *M_layer_symb=M_symb/4;
         else
             *M_layer_symb=(M_symb+2)/4;
             for(i=0;i<*M_layer_symb;i++)
             {
                 x_re[i]=d_re[4*i];
                 x_im[i]=d_im[4*i];
                 x_re[i+*M_layer_symb]=d_re[4*i+1];
                 x_im[i+*M_layer_symb]=d_im[4*i+1];
                 x_re[i+*M_layer_symb*2]=d_re[4*i+2];
                 x_im[i+*M_layer_symb*2]=d_im[4*i+2];
                 x_re[i+*M_layer_symb*3]=d_re[4*i+3];
                 x_im[i+*M_layer_symb]=d_im[4*i+3];
             }
        }


    }

    else if (N_ant==2)
    {
        if(N_codewords==1)
        {
            *M_layer_symb=M_symb/2;
            for(i=0;i<*M_layer_symb;i++)
            {
             x_re[i]=d_re[2*i];
             x_im[i]=d_im[2*i];
             x_re[i+*M_layer_symb]=d_re[2*i+1];
             x_im[i+*M_layer_symb]=d_im[2*i+1];
            }
        }
        else if(N_codewords==2)
        {
           *M_layer_symb=M_symb;
            for(i=0;i<*M_layer_symb;i++)
            {
             x_re[i]=d_re[i];
             x_im[i]=d_im[i];
            x_re[i+*M_layer_symb]=d_re[i+M_symb];
            x_im[i+*M_layer_symb]=d_im[i+M_symb];
        }
        }
    }
    else if(N_ant==3)
    {
        if(N_codewords==1)
        {
            *M_layer_symb=M_symb/3;
            for(i=0;i<*M_layer_symb;i++)
            {
             x_re[i]=d_re[3*i];
             x_im[i]=d_im[3*i];
             x_re[i+*M_layer_symb]=d_re[3*i+1];
             x_im[i+*M_layer_symb]=d_im[3*i+1];
             x_re[i+*M_layer_symb*2]=d_re[3*i+2];
             x_im[i+*M_layer_symb*2]=d_im[3*i+2];
            }
        }
        else if(N_codewords==2)
        {
           *M_layer_symb=M_symb/2;
            for(i=0;i<*M_layer_symb;i++)
            {
             x_re[i]=d_re[i];
                x_im[i]=d_im[i];
                x_re[i+*M_layer_symb]=d_re[2*i+M_symb];
                x_im[i+*M_layer_symb]=d_im[2*i+M_symb];
            x_re[i+*M_layer_symb*2]=d_re[2*i+M_symb+1];
            x_im[i+*M_layer_symb*2]=d_im[2*i+M_symb+1];
        }
        }

    }
    else if(N_ant==4)
       {
           if(N_codewords==1)
           {
               *M_layer_symb=M_symb/4;
               for(i=0;i<*M_layer_symb;i++)
               {
                x_re[i]=d_re[4*i];
                x_im[i]=d_im[4*i];
                x_re[i+*M_layer_symb]=d_re[4*i+1];
                x_im[i+*M_layer_symb]=d_im[4*i+1];
                x_re[i+*M_layer_symb*2]=d_re[4*i+2];
                x_im[i+*M_layer_symb*2]=d_im[4*i+2];
                x_re[i+*M_layer_symb*3]=d_re[4*i+3];
                x_im[i+*M_layer_symb*3]=d_im[4*i+3];
               }
           }
           else if(N_codewords==2)
           {
              *M_layer_symb=M_symb/2;
               for(i=0;i<*M_layer_symb;i++)
               {
                x_re[i]=d_re[2*i];
                x_im[i]=d_im[2*i];
               x_re[i+*M_layer_symb]=d_re[2*i+1];
               x_im[i+*M_layer_symb]=d_im[2*i+1];
               x_re[i+*M_layer_symb*2]=d_re[2*i+M_symb];
               x_im[i+*M_layer_symb*2]=d_im[2*i+M_symb];
               x_re[i+*M_layer_symb*3]=d_re[2*i+M_symb+1];
               x_im[i+*M_layer_symb*3]=d_im[2*i+M_symb+1];
               }
           }

       }
      else if(N_ant==5)
      {

              *M_layer_symb=M_symb/2;
              for(i=0;i<*M_layer_symb;i++)
              {
               x_re[i]=d_re[2*i];
               x_im[i]=d_im[2*i];
               x_re[i+*M_layer_symb]=d_re[2*i+1+M_symb];
               x_im[i+*M_layer_symb]=d_im[2*i+1+M_symb];
               x_re[i+*M_layer_symb*2]=d_re[3*i+M_symb];
               x_im[i+*M_layer_symb*2]=d_im[3*i+M_symb];
               x_re[i+*M_layer_symb*3]=d_re[3*i+1+M_symb];
               x_im[i+*M_layer_symb*3]=d_im[3*i+1+M_symb];
               x_re[i+*M_layer_symb*4]=d_re[3*i+2+M_symb];
               x_im[i+*M_layer_symb*4]=d_im[3*i+2+M_symb];
              }

      }
    else if(N_ant==6)
      {

              *M_layer_symb=M_symb/3;
              for(i=0;i<*M_layer_symb;i++)
              {
               x_re[i]=d_re[3*i];
               x_im[i]=d_im[3*i];
               x_re[i+*M_layer_symb]=d_re[3*i+1];
               x_im[i+*M_layer_symb]=d_im[3*i+1];
               x_re[i+*M_layer_symb*2]=d_re[3*i+2];
               x_im[i+*M_layer_symb*2]=d_im[3*i+2];
               x_re[i+*M_layer_symb*3]=d_re[3*i+M_symb];
               x_im[i+*M_layer_symb*3]=d_im[3*i+M_symb];
               x_re[i+*M_layer_symb*4]=d_re[3*i+1+M_symb];
               x_im[i+*M_layer_symb*4]=d_im[3*i+1+M_symb];
               x_re[i+*M_layer_symb*5]=d_re[3*i+2+M_symb];
               x_im[i+*M_layer_symb*5]=d_im[3*i+2+M_symb];
              }

      }
    else if(N_ant==7)
         {

                 *M_layer_symb=M_symb/3;
                 for(i=0;i<*M_layer_symb;i++)
                 {
                  x_re[i]=d_re[3*i];
                  x_im[i]=d_im[3*i];
                  x_re[i+*M_layer_symb]=d_re[3*i+1];
                  x_im[i+*M_layer_symb]=d_im[3*i+1];
                  x_re[i+*M_layer_symb*2]=d_re[3*i+2];
                  x_im[i+*M_layer_symb*2]=d_im[3*i+2];
                  x_re[i+*M_layer_symb*3]=d_re[4*i];
                  x_im[i+*M_layer_symb*3]=d_im[4*i];
                  x_re[i+*M_layer_symb*4]=d_re[4*i+1+M_symb];
                  x_im[i+*M_layer_symb*4]=d_im[4*i+1+M_symb];
                  x_re[i+*M_layer_symb*5]=d_re[4*i+2+M_symb];
                  x_im[i+*M_layer_symb*5]=d_im[4*i+2+M_symb];
                  x_re[i+*M_layer_symb*6]=d_re[4*i+3+M_symb];
                  x_im[i+*M_layer_symb*6]=d_im[4*i+3+M_symb];
                 }
                 }

        else if(N_ant==8)
          {

                  *M_layer_symb=M_symb/4;
                  for(i=0;i<*M_layer_symb;i++)
                  {
                   x_re[i]=d_re[4*i];
                   x_im[i]=d_im[4*i];
                   x_re[i+*M_layer_symb]=d_re[4*i+1];
                   x_im[i+*M_layer_symb]=d_im[4*i+1];
                   x_re[i+*M_layer_symb*2]=d_re[4*i+2];
                   x_im[i+*M_layer_symb*2]=d_im[4*i+2];
                   x_re[i+*M_layer_symb*3]=d_re[4*i+3];
                      x_im[i+*M_layer_symb*3]=d_im[4*i+3];
                   x_re[i+*M_layer_symb*4]=d_re[4*i+M_symb];
                   x_im[i+*M_layer_symb*4]=d_im[4*i+M_symb];
                   x_re[i+*M_layer_symb*5]=d_re[4*i+1+M_symb];
                   x_im[i+*M_layer_symb*5]=d_im[4*i+1+M_symb];
                   x_re[i+*M_layer_symb*6]=d_re[4*i+2+M_symb];
                   x_im[i+*M_layer_symb*6]=d_im[4*i+2+M_symb];
                   x_re[i+*M_layer_symb*7]=d_re[4*i+3+M_symb];
                   x_im[i+*M_layer_symb*7]=d_im[4*i+3+M_symb];
                  }
            }

    }



void modulation_mapper(int                           *bits,
                       int                           N_bits,
                       PHY_MODULATION_TYPE_ENUM  type,
                       float                           *d_re,
                       float                           *d_im,
                       int                          *M_symb)
{
    float  one_over_sqrt_2  = 1/sqrt(2);
    float  one_over_sqrt_10 = 1/sqrt(10);
    float  one_over_sqrt_42 = 1/sqrt(42);
    int i;
    int input=0;

    switch(type)
    {
    case MODULATION_TYPE_BPSK:
        // 3GPP TS 36.211 v10.1.0 section 7.1.1
        for(i=0; i<N_bits; i++)
        {
            if(0 == bits[i])
            {
                d_re[i] = one_over_sqrt_2;
                d_im[i] = one_over_sqrt_2;
            }else{
                d_re[i] = -one_over_sqrt_2;
                d_im[i] = -one_over_sqrt_2;
            }
        }
        *M_symb = N_bits;
        break;
    case MODULATION_TYPE_QPSK:
        // 3GPP TS 36.211 v10.1.0 section 7.1.2
        for(i=0; i<(N_bits/2); i++)
        {
            switch((bits[i*2] << 1) |
                   bits[i*2+1])
            {
            case 0:
                d_re[i] = +one_over_sqrt_2;
                d_im[i] = +one_over_sqrt_2;
                break;
            case 1:
                d_re[i] = +one_over_sqrt_2;
                d_im[i] = -one_over_sqrt_2;
                break;
            case 2:
                d_re[i] = -one_over_sqrt_2;
                d_im[i] = +one_over_sqrt_2;
                break;
            case 3:
                d_re[i] = -one_over_sqrt_2;
                d_im[i] = -one_over_sqrt_2;
                break;
            }
        }
        *M_symb = (N_bits/2);
        if((N_bits % 2) != 0)
        {
            *M_symb = (N_bits/2) + 1;
            // Add a trailing zero
            if(0 == bits[N_bits-1])
            {
                d_re[i] = +one_over_sqrt_2;
                d_im[i] = +one_over_sqrt_2;
            }else{
                d_re[i] = -one_over_sqrt_2;
                d_im[i] = +one_over_sqrt_2;
            }
        }
        break;
    case MODULATION_TYPE_16QAM:
        // 3GPP TS 36.211 v10.1.0 section 7.1.3
        for(i=0; i<(N_bits/4); i++)
        {
            switch((bits[i*4+0] << 3) |
                   (bits[i*4+1] << 2) |
                   (bits[i*4+2] << 1) |
                   bits[i*4+3])
            {
            case 0:
                d_re[i] = +1*one_over_sqrt_10;
                d_im[i] = +1*one_over_sqrt_10;
                break;
            case 1:
                d_re[i] = +1*one_over_sqrt_10;
                d_im[i] = +3*one_over_sqrt_10;
                break;
            case 2:
                d_re[i] = +3*one_over_sqrt_10;
                d_im[i] = +1*one_over_sqrt_10;
                break;
            case 3:
                d_re[i] = +3*one_over_sqrt_10;
                d_im[i] = +3*one_over_sqrt_10;
                break;
            case 4:
                d_re[i] = +1*one_over_sqrt_10;
                d_im[i] = -1*one_over_sqrt_10;
                break;
            case 5:
                d_re[i] = +1*one_over_sqrt_10;
                d_im[i] = -3*one_over_sqrt_10;
                break;
            case 6:
                d_re[i] = +3*one_over_sqrt_10;
                d_im[i] = -1*one_over_sqrt_10;
                break;
            case 7:
                d_re[i] = +3*one_over_sqrt_10;
                d_im[i] = -3*one_over_sqrt_10;
                break;
            case 8:
                d_re[i] = -1*one_over_sqrt_10;
                d_im[i] = +1*one_over_sqrt_10;
                break;
            case 9:
                d_re[i] = -1*one_over_sqrt_10;
                d_im[i] = +3*one_over_sqrt_10;
                break;
            case 10:
                d_re[i] = -3*one_over_sqrt_10;
                d_im[i] = +1*one_over_sqrt_10;
                break;
            case 11:
                d_re[i] = -3*one_over_sqrt_10;
                d_im[i] = +3*one_over_sqrt_10;
                break;
            case 12:
                d_re[i] = -1*one_over_sqrt_10;
                d_im[i] = -1*one_over_sqrt_10;
                break;
            case 13:
                d_re[i] = -1*one_over_sqrt_10;
                d_im[i] = -3*one_over_sqrt_10;
                break;
            case 14:
                d_re[i] = -3*one_over_sqrt_10;
                d_im[i] = -1*one_over_sqrt_10;
                break;
            case 15:
                d_re[i] = -3*one_over_sqrt_10;
                d_im[i] = -3*one_over_sqrt_10;
                break;
            }
        }
        *M_symb = (N_bits/4);
        if((N_bits % 4) != 0)
        {
            *M_symb = (N_bits/4) + 1;
            if((N_bits % 4) == 1)
            {
                input = bits[N_bits-1] << 3;
            }else if((N_bits % 4) == 2){
                input = ((bits[N_bits-2] << 3) |
                        (bits[N_bits-1] << 2));
            }else if((N_bits % 4) == 3){
                input = ((bits[N_bits-3] << 3) |
                         (bits[N_bits-2] << 2) |
                         (bits[N_bits-1] << 1));
            }
            switch(input)
            {
            case 0:
                d_re[N_bits/4] = +1*one_over_sqrt_10;
                d_im[N_bits/4] = +1*one_over_sqrt_10;
                break;
            case 1:
                d_re[N_bits/4] = +1*one_over_sqrt_10;
                d_im[N_bits/4] = +3*one_over_sqrt_10;
                break;
            case 2:
                d_re[N_bits/4] = +3*one_over_sqrt_10;
                d_im[N_bits/4] = +1*one_over_sqrt_10;
                break;
            case 3:
                d_re[N_bits/4] = +3*one_over_sqrt_10;
                d_im[N_bits/4] = +3*one_over_sqrt_10;
                break;
            case 4:
                d_re[N_bits/4] = +1*one_over_sqrt_10;
                d_im[N_bits/4] = -1*one_over_sqrt_10;
                break;
            case 5:
                d_re[N_bits/4] = +1*one_over_sqrt_10;
                d_im[N_bits/4] = -3*one_over_sqrt_10;
                break;
            case 6:
                d_re[N_bits/4] = +3*one_over_sqrt_10;
                d_im[N_bits/4] = -1*one_over_sqrt_10;
                break;
            case 7:
                d_re[N_bits/4] = +3*one_over_sqrt_10;
                d_im[N_bits/4] = -3*one_over_sqrt_10;
                break;
            case 8:
                d_re[N_bits/4] = -1*one_over_sqrt_10;
                d_im[N_bits/4] = +1*one_over_sqrt_10;
                break;
            case 9:
                d_re[N_bits/4] = -1*one_over_sqrt_10;
                d_im[N_bits/4] = +3*one_over_sqrt_10;
                break;
            case 10:
                d_re[N_bits/4] = -3*one_over_sqrt_10;
                d_im[N_bits/4] = +1*one_over_sqrt_10;
                break;
            case 11:
                d_re[N_bits/4] = -3*one_over_sqrt_10;
                d_im[N_bits/4] = +3*one_over_sqrt_10;
                break;
            case 12:
                d_re[N_bits/4] = -1*one_over_sqrt_10;
                d_im[N_bits/4] = -1*one_over_sqrt_10;
                break;
            case 13:
                d_re[N_bits/4] = -1*one_over_sqrt_10;
                d_im[N_bits/4] = -3*one_over_sqrt_10;
                break;
            case 14:
                d_re[N_bits/4] = -3*one_over_sqrt_10;
                d_im[N_bits/4] = -1*one_over_sqrt_10;
                break;
            case 15:
                d_re[N_bits/4] = -3*one_over_sqrt_10;
                d_im[N_bits/4] = -3*one_over_sqrt_10;
                break;
            }
        }
        break;
    case MODULATION_TYPE_64QAM:
        // 3GPP TS 36.211 v10.1.0 section 7.1.4
        for(i=0; i<(N_bits/6); i++)
        {
            switch((bits[i*6+0] << 5) |
                   (bits[i*6+1] << 4) |
                   (bits[i*6+2] << 3) |
                   (bits[i*6+3] << 2) |
                   (bits[i*6+4] << 1) |
                   bits[i*6+5])
            {
            case 0:
                d_re[i] = +3*one_over_sqrt_42;
                d_im[i] = +3*one_over_sqrt_42;
                break;
            case 1:
                d_re[i] = +3*one_over_sqrt_42;
                d_im[i] = +1*one_over_sqrt_42;
                break;
            case 2:
                d_re[i] = +1*one_over_sqrt_42;
                d_im[i] = +3*one_over_sqrt_42;
                break;
            case 3:
                d_re[i] = +1*one_over_sqrt_42;
                d_im[i] = +1*one_over_sqrt_42;
                break;
            case 4:
                d_re[i] = +3*one_over_sqrt_42;
                d_im[i] = +5*one_over_sqrt_42;
                break;
            case 5:
                d_re[i] = +3*one_over_sqrt_42;
                d_im[i] = +7*one_over_sqrt_42;
                break;
            case 6:
                d_re[i] = +1*one_over_sqrt_42;
                d_im[i] = +5*one_over_sqrt_42;
                break;
            case 7:
                d_re[i] = +1*one_over_sqrt_42;
                d_im[i] = +7*one_over_sqrt_42;
                break;
            case 8:
                d_re[i] = +5*one_over_sqrt_42;
                d_im[i] = +3*one_over_sqrt_42;
                break;
            case 9:
                d_re[i] = +5*one_over_sqrt_42;
                d_im[i] = +1*one_over_sqrt_42;
                break;
            case 10:
                d_re[i] = +7*one_over_sqrt_42;
                d_im[i] = +3*one_over_sqrt_42;
                break;
            case 11:
                d_re[i] = +7*one_over_sqrt_42;
                d_im[i] = +1*one_over_sqrt_42;
                break;
            case 12:
                d_re[i] = +5*one_over_sqrt_42;
                d_im[i] = +5*one_over_sqrt_42;
                break;
            case 13:
                d_re[i] = +5*one_over_sqrt_42;
                d_im[i] = +7*one_over_sqrt_42;
                break;
            case 14:
                d_re[i] = +7*one_over_sqrt_42;
                d_im[i] = +5*one_over_sqrt_42;
                break;
            case 15:
                d_re[i] = +7*one_over_sqrt_42;
                d_im[i] = +7*one_over_sqrt_42;
                break;
            case 16:
                d_re[i] = +3*one_over_sqrt_42;
                d_im[i] = -3*one_over_sqrt_42;
                break;
            case 17:
                d_re[i] = +3*one_over_sqrt_42;
                d_im[i] = -1*one_over_sqrt_42;
                break;
            case 18:
                d_re[i] = +1*one_over_sqrt_42;
                d_im[i] = -3*one_over_sqrt_42;
                break;
            case 19:
                d_re[i] = +1*one_over_sqrt_42;
                d_im[i] = -1*one_over_sqrt_42;
                break;
            case 20:
                d_re[i] = +3*one_over_sqrt_42;
                d_im[i] = -5*one_over_sqrt_42;
                break;
            case 21:
                d_re[i] = +3*one_over_sqrt_42;
                d_im[i] = -7*one_over_sqrt_42;
                break;
            case 22:
                d_re[i] = +1*one_over_sqrt_42;
                d_im[i] = -5*one_over_sqrt_42;
                break;
            case 23:
                d_re[i] = +1*one_over_sqrt_42;
                d_im[i] = -7*one_over_sqrt_42;
                break;
            case 24:
                d_re[i] = +5*one_over_sqrt_42;
                d_im[i] = -3*one_over_sqrt_42;
                break;
            case 25:
                d_re[i] = +5*one_over_sqrt_42;
                d_im[i] = -1*one_over_sqrt_42;
                break;
            case 26:
                d_re[i] = +7*one_over_sqrt_42;
                d_im[i] = -3*one_over_sqrt_42;
                break;
            case 27:
                d_re[i] = +7*one_over_sqrt_42;
                d_im[i] = -1*one_over_sqrt_42;
                break;
            case 28:
                d_re[i] = +5*one_over_sqrt_42;
                d_im[i] = -5*one_over_sqrt_42;
                break;
            case 29:
                d_re[i] = +5*one_over_sqrt_42;
                d_im[i] = -7*one_over_sqrt_42;
                break;
            case 30:
                d_re[i] = +7*one_over_sqrt_42;
                d_im[i] = -5*one_over_sqrt_42;
                break;
            case 31:
                d_re[i] = +7*one_over_sqrt_42;
                d_im[i] = -7*one_over_sqrt_42;
                break;
            case 32:
                d_re[i] = -3*one_over_sqrt_42;
                d_im[i] = +3*one_over_sqrt_42;
                break;
            case 33:
                d_re[i] = -3*one_over_sqrt_42;
                d_im[i] = +1*one_over_sqrt_42;
                break;
            case 34:
                d_re[i] = -1*one_over_sqrt_42;
                d_im[i] = +3*one_over_sqrt_42;
                break;
            case 35:
                d_re[i] = -1*one_over_sqrt_42;
                d_im[i] = +1*one_over_sqrt_42;
                break;
            case 36:
                d_re[i] = -3*one_over_sqrt_42;
                d_im[i] = +5*one_over_sqrt_42;
                break;
            case 37:
                d_re[i] = -3*one_over_sqrt_42;
                d_im[i] = +7*one_over_sqrt_42;
                break;
            case 38:
                d_re[i] = -1*one_over_sqrt_42;
                d_im[i] = +5*one_over_sqrt_42;
                break;
            case 39:
                d_re[i] = -1*one_over_sqrt_42;
                d_im[i] = +7*one_over_sqrt_42;
                break;
            case 40:
                d_re[i] = -5*one_over_sqrt_42;
                d_im[i] = +3*one_over_sqrt_42;
                break;
            case 41:
                d_re[i] = -5*one_over_sqrt_42;
                d_im[i] = +1*one_over_sqrt_42;
                break;
            case 42:
                d_re[i] = -7*one_over_sqrt_42;
                d_im[i] = +3*one_over_sqrt_42;
                break;
            case 43:
                d_re[i] = -7*one_over_sqrt_42;
                d_im[i] = +1*one_over_sqrt_42;
                break;
            case 44:
                d_re[i] = -5*one_over_sqrt_42;
                d_im[i] = +5*one_over_sqrt_42;
                break;
            case 45:
                d_re[i] = -5*one_over_sqrt_42;
                d_im[i] = +7*one_over_sqrt_42;
                break;
            case 46:
                d_re[i] = -7*one_over_sqrt_42;
                d_im[i] = +5*one_over_sqrt_42;
                break;
            case 47:
                d_re[i] = -7*one_over_sqrt_42;
                d_im[i] = +7*one_over_sqrt_42;
                break;
            case 48:
                d_re[i] = -3*one_over_sqrt_42;
                d_im[i] = -3*one_over_sqrt_42;
                break;
            case 49:
                d_re[i] = -3*one_over_sqrt_42;
                d_im[i] = -1*one_over_sqrt_42;
                break;
            case 50:
                d_re[i] = -1*one_over_sqrt_42;
                d_im[i] = -3*one_over_sqrt_42;
                break;
            case 51:
                d_re[i] = -1*one_over_sqrt_42;
                d_im[i] = -1*one_over_sqrt_42;
                break;
            case 52:
                d_re[i] = -3*one_over_sqrt_42;
                d_im[i] = -5*one_over_sqrt_42;
                break;
            case 53:
                d_re[i] = -3*one_over_sqrt_42;
                d_im[i] = -7*one_over_sqrt_42;
                break;
            case 54:
                d_re[i] = -1*one_over_sqrt_42;
                d_im[i] = -5*one_over_sqrt_42;
                break;
            case 55:
                d_re[i] = -1*one_over_sqrt_42;
                d_im[i] = -7*one_over_sqrt_42;
                break;
            case 56:
                d_re[i] = -5*one_over_sqrt_42;
                d_im[i] = -3*one_over_sqrt_42;
                break;
            case 57:
                d_re[i] = -5*one_over_sqrt_42;
                d_im[i] = -1*one_over_sqrt_42;
                break;
            case 58:
                d_re[i] = -7*one_over_sqrt_42;
                d_im[i] = -3*one_over_sqrt_42;
                break;
            case 59:
                d_re[i] = -7*one_over_sqrt_42;
                d_im[i] = -1*one_over_sqrt_42;
                break;
            case 60:
                d_re[i] = -5*one_over_sqrt_42;
                d_im[i] = -5*one_over_sqrt_42;
                break;
            case 61:
                d_re[i] = -5*one_over_sqrt_42;
                d_im[i] = -7*one_over_sqrt_42;
                break;
            case 62:
                d_re[i] = -7*one_over_sqrt_42;
                d_im[i] = -5*one_over_sqrt_42;
                break;
            case 63:
                d_re[i] = -7*one_over_sqrt_42;
                d_im[i] = -7*one_over_sqrt_42;
                break;
            }
        }
        *M_symb = (N_bits/6);
        if((N_bits % 6) != 0)
        {
            *M_symb = (N_bits/6) + 1;
            if((N_bits % 6) == 1)
            {
                input = bits[N_bits-1] << 5;
            }else if((N_bits % 6) == 2){
                input = ((bits[N_bits-2] << 5) |
                         (bits[N_bits-1] << 4));
            }else if((N_bits % 6) == 3){
                input = ((bits[N_bits-3] << 5) |
                         (bits[N_bits-2] << 4) |
                         (bits[N_bits-1] << 3));
            }else if((N_bits % 6) == 4){
                input = ((bits[N_bits-4] << 5) |
                         (bits[N_bits-3] << 4) |
                         (bits[N_bits-2] << 3) |
                         (bits[N_bits-1] << 2));
            }else if((N_bits % 6) == 5){
                input = ((bits[N_bits-5] << 5) |
                         (bits[N_bits-4] << 4) |
                         (bits[N_bits-3] << 3) |
                         (bits[N_bits-2] << 2) |
                         (bits[N_bits-1] << 1));
            }
            switch(input)
            {
            case 0:
                d_re[N_bits/6] = +3*one_over_sqrt_42;
                d_im[N_bits/6] = +3*one_over_sqrt_42;
                break;
            case 1:
                d_re[N_bits/6] = +3*one_over_sqrt_42;
                d_im[N_bits/6] = +1*one_over_sqrt_42;
                break;
            case 2:
                d_re[N_bits/6] = +1*one_over_sqrt_42;
                d_im[N_bits/6] = +3*one_over_sqrt_42;
                break;
            case 3:
                d_re[N_bits/6] = +1*one_over_sqrt_42;
                d_im[N_bits/6] = +1*one_over_sqrt_42;
                break;
            case 4:
                d_re[N_bits/6] = +3*one_over_sqrt_42;
                d_im[N_bits/6] = +5*one_over_sqrt_42;
                break;
            case 5:
                d_re[N_bits/6] = +3*one_over_sqrt_42;
                d_im[N_bits/6] = +7*one_over_sqrt_42;
                break;
            case 6:
                d_re[N_bits/6] = +1*one_over_sqrt_42;
                d_im[N_bits/6] = +5*one_over_sqrt_42;
                break;
            case 7:
                d_re[N_bits/6] = +1*one_over_sqrt_42;
                d_im[N_bits/6] = +7*one_over_sqrt_42;
                break;
            case 8:
                d_re[N_bits/6] = +5*one_over_sqrt_42;
                d_im[N_bits/6] = +3*one_over_sqrt_42;
                break;
            case 9:
                d_re[N_bits/6] = +5*one_over_sqrt_42;
                d_im[N_bits/6] = +1*one_over_sqrt_42;
                break;
            case 10:
                d_re[N_bits/6] = +7*one_over_sqrt_42;
                d_im[N_bits/6] = +3*one_over_sqrt_42;
                break;
            case 11:
                d_re[N_bits/6] = +7*one_over_sqrt_42;
                d_im[N_bits/6] = +1*one_over_sqrt_42;
                break;
            case 12:
                d_re[N_bits/6] = +5*one_over_sqrt_42;
                d_im[N_bits/6] = +5*one_over_sqrt_42;
                break;
            case 13:
                d_re[N_bits/6] = +5*one_over_sqrt_42;
                d_im[N_bits/6] = +7*one_over_sqrt_42;
                break;
            case 14:
                d_re[N_bits/6] = +7*one_over_sqrt_42;
                d_im[N_bits/6] = +5*one_over_sqrt_42;
                break;
            case 15:
                d_re[N_bits/6] = +7*one_over_sqrt_42;
                d_im[N_bits/6] = +7*one_over_sqrt_42;
                break;
            case 16:
                d_re[N_bits/6] = +3*one_over_sqrt_42;
                d_im[N_bits/6] = -3*one_over_sqrt_42;
                break;
            case 17:
                d_re[N_bits/6] = +3*one_over_sqrt_42;
                d_im[N_bits/6] = -1*one_over_sqrt_42;
                break;
            case 18:
                d_re[N_bits/6] = +1*one_over_sqrt_42;
                d_im[N_bits/6] = -3*one_over_sqrt_42;
                break;
            case 19:
                d_re[N_bits/6] = +1*one_over_sqrt_42;
                d_im[N_bits/6] = -1*one_over_sqrt_42;
                break;
            case 20:
                d_re[N_bits/6] = +3*one_over_sqrt_42;
                d_im[N_bits/6] = -5*one_over_sqrt_42;
                break;
            case 21:
                d_re[N_bits/6] = +3*one_over_sqrt_42;
                d_im[N_bits/6] = -7*one_over_sqrt_42;
                break;
            case 22:
                d_re[N_bits/6] = +1*one_over_sqrt_42;
                d_im[N_bits/6] = -5*one_over_sqrt_42;
                break;
            case 23:
                d_re[N_bits/6] = +1*one_over_sqrt_42;
                d_im[N_bits/6] = -7*one_over_sqrt_42;
                break;
            case 24:
                d_re[N_bits/6] = +5*one_over_sqrt_42;
                d_im[N_bits/6] = -3*one_over_sqrt_42;
                break;
            case 25:
                d_re[N_bits/6] = +5*one_over_sqrt_42;
                d_im[N_bits/6] = -1*one_over_sqrt_42;
                break;
            case 26:
                d_re[N_bits/6] = +7*one_over_sqrt_42;
                d_im[N_bits/6] = -3*one_over_sqrt_42;
                break;
            case 27:
                d_re[N_bits/6] = +7*one_over_sqrt_42;
                d_im[N_bits/6] = -1*one_over_sqrt_42;
                break;
            case 28:
                d_re[N_bits/6] = +5*one_over_sqrt_42;
                d_im[N_bits/6] = -5*one_over_sqrt_42;
                break;
            case 29:
                d_re[N_bits/6] = +5*one_over_sqrt_42;
                d_im[N_bits/6] = -7*one_over_sqrt_42;
                break;
            case 30:
                d_re[N_bits/6] = +7*one_over_sqrt_42;
                d_im[N_bits/6] = -5*one_over_sqrt_42;
                break;
            case 31:
                d_re[N_bits/6] = +7*one_over_sqrt_42;
                d_im[N_bits/6] = -7*one_over_sqrt_42;
                break;
            case 32:
                d_re[N_bits/6] = -3*one_over_sqrt_42;
                d_im[N_bits/6] = +3*one_over_sqrt_42;
                break;
            case 33:
                d_re[N_bits/6] = -3*one_over_sqrt_42;
                d_im[N_bits/6] = +1*one_over_sqrt_42;
                break;
            case 34:
                d_re[N_bits/6] = -1*one_over_sqrt_42;
                d_im[N_bits/6] = +3*one_over_sqrt_42;
                break;
            case 35:
                d_re[N_bits/6] = -1*one_over_sqrt_42;
                d_im[N_bits/6] = +1*one_over_sqrt_42;
                break;
            case 36:
                d_re[N_bits/6] = -3*one_over_sqrt_42;
                d_im[N_bits/6] = +5*one_over_sqrt_42;
                break;
            case 37:
                d_re[N_bits/6] = -3*one_over_sqrt_42;
                d_im[N_bits/6] = +7*one_over_sqrt_42;
                break;
            case 38:
                d_re[N_bits/6] = -1*one_over_sqrt_42;
                d_im[N_bits/6] = +5*one_over_sqrt_42;
                break;
            case 39:
                d_re[N_bits/6] = -1*one_over_sqrt_42;
                d_im[N_bits/6] = +7*one_over_sqrt_42;
                break;
            case 40:
                d_re[N_bits/6] = -5*one_over_sqrt_42;
                d_im[N_bits/6] = +3*one_over_sqrt_42;
                break;
            case 41:
                d_re[N_bits/6] = -5*one_over_sqrt_42;
                d_im[N_bits/6] = +1*one_over_sqrt_42;
                break;
            case 42:
                d_re[N_bits/6] = -7*one_over_sqrt_42;
                d_im[N_bits/6] = +3*one_over_sqrt_42;
                break;
            case 43:
                d_re[N_bits/6] = -7*one_over_sqrt_42;
                d_im[N_bits/6] = +1*one_over_sqrt_42;
                break;
            case 44:
                d_re[N_bits/6] = -5*one_over_sqrt_42;
                d_im[N_bits/6] = +5*one_over_sqrt_42;
                break;
            case 45:
                d_re[N_bits/6] = -5*one_over_sqrt_42;
                d_im[N_bits/6] = +7*one_over_sqrt_42;
                break;
            case 46:
                d_re[N_bits/6] = -7*one_over_sqrt_42;
                d_im[N_bits/6] = +5*one_over_sqrt_42;
                break;
            case 47:
                d_re[N_bits/6] = -7*one_over_sqrt_42;
                d_im[N_bits/6] = +7*one_over_sqrt_42;
                break;
            case 48:
                d_re[N_bits/6] = -3*one_over_sqrt_42;
                d_im[N_bits/6] = -3*one_over_sqrt_42;
                break;
            case 49:
                d_re[N_bits/6] = -3*one_over_sqrt_42;
                d_im[N_bits/6] = -1*one_over_sqrt_42;
                break;
            case 50:
                d_re[N_bits/6] = -1*one_over_sqrt_42;
                d_im[N_bits/6] = -3*one_over_sqrt_42;
                break;
            case 51:
                d_re[N_bits/6] = -1*one_over_sqrt_42;
                d_im[N_bits/6] = -1*one_over_sqrt_42;
                break;
            case 52:
                d_re[N_bits/6] = -3*one_over_sqrt_42;
                d_im[N_bits/6] = -5*one_over_sqrt_42;
                break;
            case 53:
                d_re[N_bits/6] = -3*one_over_sqrt_42;
                d_im[N_bits/6] = -7*one_over_sqrt_42;
                break;
            case 54:
                d_re[N_bits/6] = -1*one_over_sqrt_42;
                d_im[N_bits/6] = -5*one_over_sqrt_42;
                break;
            case 55:
                d_re[N_bits/6] = -1*one_over_sqrt_42;
                d_im[N_bits/6] = -7*one_over_sqrt_42;
                break;
            case 56:
                d_re[N_bits/6] = -5*one_over_sqrt_42;
                d_im[N_bits/6] = -3*one_over_sqrt_42;
                break;
            case 57:
                d_re[N_bits/6] = -5*one_over_sqrt_42;
                d_im[N_bits/6] = -1*one_over_sqrt_42;
                break;
            case 58:
                d_re[N_bits/6] = -7*one_over_sqrt_42;
                d_im[N_bits/6] = -3*one_over_sqrt_42;
                break;
            case 59:
                d_re[N_bits/6] = -7*one_over_sqrt_42;
                d_im[N_bits/6] = -1*one_over_sqrt_42;
                break;
            case 60:
                d_re[N_bits/6] = -5*one_over_sqrt_42;
                d_im[N_bits/6] = -5*one_over_sqrt_42;
                break;
            case 61:
                d_re[N_bits/6] = -5*one_over_sqrt_42;
                d_im[N_bits/6] = -7*one_over_sqrt_42;
                break;
            case 62:
                d_re[N_bits/6] = -7*one_over_sqrt_42;
                d_im[N_bits/6] = -5*one_over_sqrt_42;
                break;
            case 63:
                d_re[N_bits/6] = -7*one_over_sqrt_42;
                d_im[N_bits/6] = -7*one_over_sqrt_42;
                break;
            }
        }
        break;
    }
}

void calc_crc(int  *a_bits,
              int  N_a_bits,
              int  crc,
              int  *p_bits,
              int  N_p_bits)
{
    int i;
    int crc_rem   = 0;
    int crc_check = (1 << N_p_bits);
    int  tmp_array[N_a_bits + N_p_bits];
    tmp_array[N_a_bits]= a_bits[N_a_bits];

    // Initialize tmp_array


    for(i=0; i<N_a_bits + N_p_bits; i++)
    {
        crc_rem <<= 1;
        crc_rem  |= tmp_array[i];

        if(crc_rem & crc_check)
            crc_rem ^= crc;

    }

    for(i=0; i<N_p_bits; i++)
    {
        p_bits[i] = (crc_rem >> (N_p_bits-1-i)) & 1;

    }
}


void rate_match_conv(LIBLTE_PHY *phy_struct,
int             *d_bits,
int             N_d_bits,
int             N_e_bits,
int             *e_bits)
{
    int R_c;
    int C_c=32;
    int i,j,k,x,N_dummy;
    int i_dx;
    int w_idx=0;
    int K_pi,K_w;
    R_c=0;

    while((N_d_bits/3)>R_c*C_c)
    {
        R_c++;
    }

    for(x=0;x<3;x++)
    {
    if(R_c*C_c>(N_d_bits/3))
        N_dummy=R_c*C_c-N_d_bits/3;
    else
        N_dummy=0;
    for(i=0;i<N_dummy;i++)
        phy_struct->rmc_tmp[i]=100;
    i_dx=0;
    for(i=N_dummy;i<C_c*R_c;i++)
    {
        phy_struct->rmc_tmp[i]=d_bits[i_dx*3+x];
        i_dx++;
    }



        i_dx=0;
    for(i=0;i<R_c;i++)
    { for(j=0;j<C_c;j++)
      {
        phy_struct->rmc_sb_mat[i][j]=phy_struct->rmc_tmp[i_dx++];


       }

    }

    //Permutation
    for(i=0;i<R_c;i++)
    {
        for(j=0;j<C_c;j++)
        {
            phy_struct->rmc_sb_perm_mat[i][j]=phy_struct->rmc_sb_mat[i][ICC_PERM[j]];

        }

    }
    for(j=0;j<C_c;j++)
    {
        for(i=0;i<R_c;i++)
        {
        phy_struct->rmc_w[w_idx++]=phy_struct->rmc_sb_perm_mat[i][j];
        }


        }


}
 //Creating the circular buffer
    K_pi = R_c*C_c;
    K_w=3*K_pi;
    k=0;
    j=0;
    while(k<N_e_bits)
    {
        if(phy_struct->rmc_w[j%K_w]!=100)
        {
            e_bits[k++]=phy_struct->rmc_w[j%K_w];

        }
        j++;
    }

}


void generate_crs(int ns,int l,bool x,int N_id_cell,float *c_re,float *c_im)
{
    int c_init,i;
    int n_sbar,N_cp;
    int len= 2*N_rb_max;
    int c[2*len];
    if(x==0)
    {
        n_sbar=ns;
        N_cp=1;
    }
    else
    {
        n_sbar=10*ceil(ns/10)+ns%2;
        N_cp=0;
    }

    c_init=((7*(n_sbar+1)+l+1)*(2*N_id_cell+1)>>10 )+2*N_id_cell+N_cp;

    scramb_generate(c_init,len,c);

    for(i=0;i<len;i++)
    {
        c_re[i]=(1/sqrt(2))*(1-2*c[2*i]);
        c_im[i]=(1/sqrt(2))*(1-2*c[2*i+1]);
    }

}
void pbch_channel_map(FRAME *subframe,LIBLTE_PHY *phy_struct,int *in_bits,int sfn,int N_ant,int N_id_cell)
{
       int ant_mask1[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
       int ant_mask2[16]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
       int ant_mask3[16]={0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1};
       int p,k;
       int *ant_mask=NULL;
       int *a_bits;
       int p_bits[16];
       int i;
       int N_d_bits;
       int offset;
       int M_symb;
       int M_layer_symb;
       int M_ap_symb;
       int idx;
       //int out_bits[1920];
       int g[3] = {0133, 0171, 0165};

       if (N_ant==1)
           ant_mask=ant_mask1;
       else if (N_ant==2)
           ant_mask=ant_mask2;
       else
           ant_mask=ant_mask3;

        a_bits=in_bits;
      //CRC calculation and assigning the bits
       calc_crc(&a_bits[0],24,CRC16,&p_bits[0],16);
       for(i=0;i<16;i++)
           p_bits[i]=p_bits[i]^ant_mask[i];
       for(i=0;i<24;i++)
           phy_struct->bch_bits[i]=a_bits[i];
       for(i=0;i<16;i++)
           phy_struct->bch_bits[i+24]=p_bits[i];

    if(sfn%4==0)  //Limiting the encoding to only when sfn%4=0
    {
       conv_encode(phy_struct,&phy_struct->bch_bits[0],40,7,3,g,1,&phy_struct->bch_tx_dbits[0],&N_d_bits);

       rate_match_conv(phy_struct,&phy_struct->bch_tx_dbits[0],N_d_bits,1920,phy_struct->bch_encode_bits);
       phy_struct->bch_N_bits=1920;
       scramb_generate(N_id_cell,phy_struct->bch_N_bits,&phy_struct->bch_c[0]);
    }

    //for(i=0;i<1920;i++)
    //    cout<<phy_struct->bch_encode_bits[i];
    //cout<<endl;
    offset=(sfn%4)*480;
     //cout<<"offset :"<<offset<<endl;
     //for(i=0;i<480;i++)
     //{
     //    cout<<phy_struct->bch_encode_bits[i+offset];
    //}
    //cout<<endl;


       for(i=0;i<480;i++)
       {
           phy_struct->bch_scramb_bits[i]=phy_struct->bch_c[i+offset]^phy_struct->bch_encode_bits[i+offset];


       }
       modulation_mapper(&phy_struct->bch_scramb_bits[0],480,MODULATION_TYPE_QPSK,&phy_struct->bch_d_re[0],&phy_struct->bch_d_im[0],&M_symb);


       layer_mapper(&phy_struct->bch_d_re[0],&phy_struct->bch_d_im[0],M_symb,N_ant,1,transmit,&phy_struct->bch_x_re[0],&phy_struct->bch_x_im[0],&M_layer_symb);


       pre_coder_dl(&phy_struct->bch_d_re[0],&phy_struct->bch_d_im[0],M_layer_symb,N_ant,transmit,phy_struct->bch_y_re[0] ,phy_struct->bch_y_im[0],240,&M_ap_symb);




       for(p=0;p<N_ant;p++)
       { idx=0;
           for(i=0;i<phy_struct->N_sc_rb*phy_struct->N_rb_dl;i++)
           {   k=(phy_struct->N_sc_rb*phy_struct->N_rb_dl)/2-36+i;
               if(i%3!=(N_id_cell)%3)
               {
                   subframe->tx_symb_re[p][7][k]=phy_struct->bch_y_re[p][idx];
                   subframe->tx_symb_im[p][7][k]=phy_struct->bch_y_im[p][idx];
                   subframe->tx_symb_re[p][8][k]=phy_struct->bch_y_re[p][idx+48];
                   subframe->tx_symb_im[p][8][k]=phy_struct->bch_y_re[p][idx+48];
                   idx++;
               }
               subframe->tx_symb_re[p][9][k]=phy_struct->bch_y_re[p][i+96];
               subframe->tx_symb_im[p][9][k]=phy_struct->bch_y_re[p][i+96];
               subframe->tx_symb_re[p][10][k]=phy_struct->bch_y_re[p][i+168];
               subframe->tx_symb_im[p][10][k]=phy_struct->bch_y_re[p][i+168];

               }


       }



}


void pcfich_channel_map(LIBLTE_PHY *phy_struct,int N_id_cell,int N_ant,PHY_PCFICH *pcfich,FRAME *subframe)
{
       int M_symb,M_layer_symb,M_ap_symb;
       int c_init,i;
       int l_prime;
       int p,idx,j;
       int N_bits,k_bar;

          l_prime=0;

           cfi_encode(phy_struct,pcfich->cfi,phy_struct->pdcch_encode_bits,&N_bits);
        c_init=((subframe->num)+1)*((2*N_id_cell+1)<<9) +N_id_cell;


           scramb_generate(c_init,N_bits,phy_struct->pdcch_c);

           for(i=0;i<N_bits;i++)
               phy_struct->pdcch_scramb_bits[i]=phy_struct->pdcch_encode_bits[i]^phy_struct->pdcch_c[i];


           modulation_mapper(phy_struct->pdcch_scramb_bits,32,MODULATION_TYPE_QPSK,phy_struct->pdcch_d_re,phy_struct->pdcch_d_im,&M_symb);


           layer_mapper(phy_struct->pdcch_d_re,phy_struct->pdcch_d_im,M_symb,N_ant,1,transmit,phy_struct->pdcch_x_re,phy_struct->pdcch_x_im,&M_layer_symb);




           pre_coder_dl(phy_struct->pdcch_x_re,phy_struct->pdcch_x_im,M_layer_symb,N_ant,transmit,phy_struct->pdcch_y_re[0],phy_struct->pdcch_y_im[0],576,&M_ap_symb);




       k_bar=(phy_struct->N_sc_rb/2)*(N_id_cell%(2*phy_struct->N_rb_dl));


       for(i=0;i<4;i++)
       {
           pcfich->k[i]=(k_bar+(i*phy_struct->N_rb_dl/2)*(phy_struct->N_sc_rb/2))%((phy_struct->N_sc_rb)*(phy_struct->N_rb_dl));

           pcfich->n[i]=(pcfich->k[i]/6)-0.5;



           for(p=0;p<N_ant;p++)
           {
               idx=0;
               for(j=0;j<6;j++)
               {
                   if((N_id_cell%3)!=(j%3))
                   {
                       subframe->tx_symb_re[p][0][pcfich->k[i]+j]=phy_struct->pdcch_y_re[p][4*i+idx];
                       subframe->tx_symb_im[p][0][pcfich->k[i]+j]=phy_struct->pdcch_y_im[p][4*i+idx];

                       idx++;
                   }
               }
           }
       }

}

    void phich_channel_map(LIBLTE_PHY *phy_struct,int N_id_cell,int N_ant,PHY_PHICH *PHICH,PHY_PCFICH *pcfich,FRAME *subframe)
    {
        int seq,M_symb,M_layer_symb,M_ap_symb;
           int w_idx,z_idx;
           int *vseq;
           int ack_seq[3]={1,1,1};
           int nack_seq[3]={0,0,0};
           int m_0;
           int c_init,i;
           int l_prime,n_l_prime;
           int n_hat[3];
           int y_idx,p,idx,j;


              l_prime=0;

        //No of REG's
        PHICH->N_reg=phy_struct->N_group_phich*3;

        //To generate the scrambling sequence
            c_init=(subframe->num)+1*((2*N_id_cell+1)<<9) +N_id_cell;
            scramb_generate(c_init,12,phy_struct->pdcch_c);

            idx=0;
          //To generate the 12 bit sequence and multiplexing
            for(m_0=0;m_0<phy_struct->N_group_phich;m_0++)
            {
                for(i=0;i<12;i++)
                {
                    phy_struct->pdcch_d_re[i]=0;
                    phy_struct->pdcch_d_im[i]=0;

                    }
                for(seq=0;seq<8;seq++)
                {
                    if(PHICH->present[m_0][seq])
                     {if(PHICH->b[m_0][seq])
                        vseq =ack_seq;
                        else
                            vseq =nack_seq;

                    modulation_mapper(&vseq[0],3,MODULATION_TYPE_BPSK,&PHICH->z_re[0],&PHICH->z_im[0],&M_symb);

                    for(i=0;i<phy_struct->N_sf_phich*3;i++)
                    {
                        w_idx= i%phy_struct->N_sf_phich;
                        z_idx=i/phy_struct->N_sf_phich;



                        if(phy_struct->pdcch_c[i])
                           {
                            phy_struct->pdcch_d_re[i]=PHICH_normal_re[seq][w_idx]*(-PHICH->z_re[z_idx])+PHICH_normal_im[seq][w_idx]*(PHICH->z_im[z_idx]);
                        phy_struct->pdcch_d_im[i]=PHICH_normal_re[seq][w_idx]*(-PHICH->z_im[z_idx])+PHICH_normal_im[seq][w_idx]*(-PHICH->z_re[z_idx]);
                            }

                        else{ phy_struct->pdcch_d_re[i]=PHICH_normal_re[seq][w_idx]*(PHICH->z_re[z_idx])+PHICH_normal_im[seq][w_idx]*(-PHICH->z_im[z_idx]);
                        phy_struct->pdcch_d_im[i]=PHICH_normal_re[seq][w_idx]*(PHICH->z_im[z_idx])+PHICH_normal_im[seq][w_idx]*(PHICH->z_re[z_idx]);
                        }


                    }
                    }
                }




            //Layer mapping

            layer_mapper(phy_struct->pdcch_d_re,phy_struct->pdcch_d_im,3*phy_struct->N_sf_phich,N_ant,1,transmit,&phy_struct->pdcch_x_re[0],&phy_struct->pdcch_x_im[0],&M_layer_symb);

            //Precoding

             pre_coder_dl(&phy_struct->pdcch_x_re[0],&phy_struct->pdcch_x_im[0],M_layer_symb,N_ant,transmit,phy_struct->pdcch_y_re[0] ,phy_struct->pdcch_y_im[0],576,&M_ap_symb);




            //Mapping


            n_l_prime=phy_struct->N_rb_dl*2 -4;

            for(i=0;i<3;i++)
            {
                n_hat[i]=(N_id_cell+m_0+(i*n_l_prime)/3)%n_l_prime;


            }
                for(i=0;i<pcfich->N_reg;i++)
                { for(j=0;j<3;j++)
                   {
                    if(n_hat[j]>pcfich->n[i])
                        n_hat[j]++;

                    }
                }


                for(j=0;j<4;j++)
                {for(i=0;i<3;i++)
                {
                    if(n_hat[i]*6==pcfich->k[j])
                        n_hat[i]+=1;
                }
                    }
               //to be checked

            y_idx=0;
            for(i=0;i<3;i++)
            {
                PHICH->k[idx+i]=n_hat[i]*6;

                for(p=0;p<N_ant;p++)
                {  y_idx=4*i;
                    for(j=0;j<6;j++)
                    {
                        if(N_id_cell%3 !=j%3)
                        {
                            subframe->tx_symb_re[p][0][PHICH->k[idx+i]+j]=phy_struct->pdcch_y_re[p][y_idx];
                            subframe->tx_symb_im[p][0][PHICH->k[idx+i]+j]=phy_struct->pdcch_y_im[p][y_idx];
                            y_idx++;
                        }
                    }
                }

            }
            idx+=3;

        }
        }


 void init_phy(LIBLTE_PHY *phy_struct,FS_enum fs,int N_ant,int N_id_cell,float phich_res,bool x)
{    int i;
    if(phy_struct!=NULL)
    {
    switch(fs)
    {
        case 0:
        phy_struct-> fs=1920000;
        phy_struct-> N_samps_per_symb=N_samps_per_symbol_1_92;
        phy_struct-> N_samps_cp_1_else=N_samps_CP_else_1_92;
        phy_struct-> N_samps_cp_1_0=N_samps_CP_1_92;
        phy_struct->N_samps_per_slot=N_samples_slot_1_92;
        phy_struct->N_samps_per_subfr=N_samps_per_subframe_1_92;
        phy_struct->N_samps_per_frame=N_samps_per_frame_1_92;
            phy_struct->N_rb_dl=6;
            break;
        case 1:
        phy_struct-> fs=3840000;
        phy_struct-> N_samps_per_symb=N_samps_per_symbol_3_84;
        phy_struct-> N_samps_cp_1_else=N_samps_CP_else_3_84;
        phy_struct-> N_samps_cp_1_0=N_samps_CP_3_84;
        phy_struct->N_samps_per_slot=N_samples_slot_3_84;
        phy_struct->N_samps_per_subfr=N_samps_per_subframe_3_84;
        phy_struct->N_samps_per_frame=N_samps_per_frame_3_84;
            phy_struct->N_rb_dl=15;
            break;
        case 2:
            phy_struct-> fs=7680000;
            phy_struct-> N_samps_per_symb=N_samps_per_symbol_7_68;
            phy_struct-> N_samps_cp_1_else=N_samps_CP_else_7_68;
            phy_struct-> N_samps_cp_1_0=N_samps_CP_7_68;
            phy_struct->N_samps_per_slot=N_samples_slot_7_68;
            phy_struct->N_samps_per_subfr=N_samps_per_subframe_7_68;
            phy_struct->N_samps_per_frame=N_samps_per_frame_7_68;
            phy_struct->N_rb_dl=25;
            break;
        case 3:
            phy_struct-> fs=15360000;
            phy_struct-> N_samps_per_symb=N_samps_per_symbol_15_36;
            phy_struct-> N_samps_cp_1_else=N_samps_CP_else_15_36;
            phy_struct-> N_samps_cp_1_0=N_samps_CP_15_36;
            phy_struct->N_samps_per_slot=N_samples_slot_15_36;
            phy_struct->N_samps_per_subfr=N_samps_per_subframe_15_36;
            phy_struct->N_samps_per_frame=N_samps_per_frame_15_36;
             phy_struct->N_rb_dl=50;
            break;
        case 4:
            phy_struct-> fs=30720000;
            phy_struct-> N_samps_per_symb=N_samps_per_symbol_30_72;
            phy_struct-> N_samps_cp_1_else=N_samps_CP_else_30_72;
            phy_struct-> N_samps_cp_1_0=N_samps_CP_30_72;
            phy_struct->N_samps_per_slot=N_samples_slot_30_72;
            phy_struct->N_samps_per_subfr=N_samps_per_subframe_30_72;
            phy_struct->N_samps_per_frame=N_samps_per_frame_30_72;
             phy_struct->N_rb_dl=75;
            break;
      }
        phy_struct->N_ant=N_ant;
        phy_struct->N_id_cell=N_id_cell;


        //PHICH config

        //Normal CP
        if(x==0)
        {
        phy_struct->N_group_phich=int(ceil((phich_res*phy_struct->N_rb_dl/8)));
        phy_struct->N_sf_phich=4;
        }
        //Extended CP
        else
        {
            phy_struct->N_group_phich=int(ceil((2*phich_res*phy_struct->N_rb_dl/8)));
            phy_struct->N_sf_phich=2;
        }

        //CRS generation
        phy_struct->N_id_cell_crs=N_id_cell;
        for(i=0;i<20;i++)
        {
            generate_crs(i,0,x,N_id_cell,phy_struct->crs_re_storage[i][0],phy_struct->crs_im_storage[i][0]);
            generate_crs(i,4,x,N_id_cell,phy_struct->crs_re_storage[i][1],phy_struct->crs_im_storage[i][1]);
            generate_crs(i,7,x,N_id_cell,phy_struct->crs_re_storage[i][2],phy_struct->crs_im_storage[i][2]);

        }
    }



}
void crs_map(LIBLTE_PHY *phy_struct,FRAME *subframe,int N_id_cell,int N_ant,bool x)
{
    int p,N_sym,i,j,k,m_prime;
    int sym[4];
    int v_shift=N_id_cell%6;
    int v[4];
    float *crs_re[14];
    float *crs_im[14];
    if(phy_struct->N_id_cell_crs==N_id_cell)
    {
        crs_re[0]=&phy_struct->crs_re_storage[subframe->num*2][0][0];
        crs_im[0]=&phy_struct->crs_im_storage[subframe->num*2][0][0];
        crs_re[1]=&phy_struct->crs_re_storage[subframe->num*2][1][0];
        crs_im[1]=&phy_struct->crs_im_storage[subframe->num*2][1][0];
        crs_re[4]=&phy_struct->crs_re_storage[subframe->num*2][2][0];
        crs_im[4]=&phy_struct->crs_im_storage[subframe->num*2][2][0];
        crs_re[7]=&phy_struct->crs_re_storage[subframe->num*2+1][0][0];
        crs_im[7]=&phy_struct->crs_im_storage[subframe->num*2+1][0][0];
        crs_re[8]=&phy_struct->crs_re_storage[subframe->num*2+1][1][0];
        crs_im[8]=&phy_struct->crs_im_storage[subframe->num*2+1][1][0];
        crs_re[11]=&phy_struct->crs_re_storage[subframe->num*2+1][2][0];
        crs_im[11]=&phy_struct->crs_im_storage[subframe->num*2+1][2][0];
    }
    else
    {
        generate_crs(subframe->num*2,0,x,N_id_cell,phy_struct->crs_re[0],phy_struct->crs_im[0]);
        generate_crs(subframe->num*2,1,x,N_id_cell,phy_struct->crs_re[1],phy_struct->crs_im[1]);
        generate_crs(subframe->num*2,4,x,N_id_cell,phy_struct->crs_re[4],phy_struct->crs_im[4]);
        generate_crs(subframe->num*2,7,x,N_id_cell,phy_struct->crs_re[7],phy_struct->crs_im[7]);
        generate_crs(subframe->num*2,8,x,N_id_cell,phy_struct->crs_re[8],phy_struct->crs_im[8]);
        generate_crs(subframe->num*2,11,x,N_id_cell,phy_struct->crs_re[11],phy_struct->crs_im[11]);

        crs_re[0]=&phy_struct->crs_re[0][0];
        crs_im[0]=&phy_struct->crs_im[0][0];
        crs_re[1]=&phy_struct->crs_re[1][0];
        crs_im[1]=&phy_struct->crs_im[1][0];
        crs_re[4]=&phy_struct->crs_re[4][0];
        crs_re[4]=&phy_struct->crs_im[4][0];
        crs_re[7]=&phy_struct->crs_re[7][0];
        crs_im[7]=&phy_struct->crs_im[7][0];
        crs_re[8]=&phy_struct->crs_re[8][0];
        crs_im[8]=&phy_struct->crs_im[8][0];
        crs_re[11]=&phy_struct->crs_re[11][0];
        crs_im[11]=&phy_struct->crs_im[11][0];

    }

    for(p=0;p<N_ant;p++)
    {
        if(p==0)
        {
            v[0]=0;
            v[1]=3;
            v[2]=0;
            v[3]=3;
            sym[0]=0;
            sym[1]=4;
            sym[2]=7;
            sym[3]=11;
            N_sym=4;

        }
        else if(p==1)
        {
            v[0]=3;
            v[1]=0;
            v[2]=3;
            v[3]=0;
            sym[0]=0;
            sym[1]=4;
            sym[2]=7;
            sym[3]=11;
            N_sym=4;

        }
        else if(p==2)
        {
            v[0]=0;
            v[1]=3;
            sym[0]=1;
            sym[1]=8;
            N_sym=2;

        }
        else if(p==3)
        {
            v[0]=3;
            v[1]=6;
            sym[0]=1;
            sym[1]=8;
            N_sym=2;

        }

    for(i=0;i<N_sym;i++)
    {
        for(j=0;j<2*phy_struct->N_rb_dl;j++)
        {
            k=6*j+(v[i]+v_shift)%6;

            m_prime=j+N_rb_max-phy_struct->N_rb_dl;

            subframe->tx_symb_re[p][sym[i]][k]=crs_re[sym[i]][m_prime];

            subframe->tx_symb_im[p][sym[i]][k]=crs_im[sym[i]][m_prime];



        }

    }
}
}
void dci_1a_pack(PHY_ALLOCATION_STRUCT *alloc,int N_rb_dl,int N_ant,int *out_bits,int *N_out_bits)

{   int RIV_length,RIV,N_prb_1a,size;
    int *dci=out_bits;
    //CA is not present
    convert_to_bits(0,&dci,3);
    //DCI 1A flag
    convert_to_bits( DCI_0_1A_FLAG_1A,&dci,1);
    //Check if it is a PDCCH random access procedure
    if(MAC_PRNTI==alloc->rnti||MAC_SIRNTI==alloc->rnti||(MAC_SI_RNTI_START<=alloc->rnti&&MAC_SI_RNTI_END>=alloc->rnti))
    {
        convert_to_bits(DCI_VRB_TYPE_LOCALIZED,&dci,1); //Only supporting localised should make it Distributed

        RIV_length=ceil(logf(N_rb_dl*(N_rb_dl+1))/logf(2));

        if((alloc->N_prb-1)<N_rb_dl/2)
        {
            RIV=N_rb_dl*(alloc->N_prb-1)+alloc->prb[0][0];
        }
        else
        {
           RIV=N_rb_dl*(alloc->N_prb+1)+(N_rb_dl-1-alloc->prb[0][0]);

         }

        convert_to_bits(RIV,&dci,RIV_length);
        //Modulation scheme
        convert_to_bits(alloc->mcs,&dci,5);
        //HARQ indicator
        convert_to_bits(0,&dci,3);//Only considering FDD
        //New data indicator
        convert_to_bits(0,&dci,1);
        //Redundancy version
        convert_to_bits(alloc->rv_idx,&dci,2);
        //TPC
        N_prb_1a=3;

        convert_to_bits(1,&dci,2);

        alloc->tbs=TBS[alloc->mcs][N_prb_1a-1];  //To be checked
    }
    else
    {
            convert_to_bits(DCI_VRB_TYPE_LOCALIZED,&dci,1);

            RIV_length=ceil(logf(N_rb_dl*(N_rb_dl+1))/logf(2));

            if((alloc->N_prb-1)<N_rb_dl/2)
            {
                RIV=N_rb_dl*(alloc->N_prb-1)+alloc->prb[0][0];
            }
            else
            {
               RIV=N_rb_dl*(alloc->N_prb+1)+(N_rb_dl-1-alloc->prb[0][0]);

             }

            convert_to_bits(RIV,&dci,RIV_length);
            //Modulation scheme
            convert_to_bits(alloc->mcs,&dci,5);
            //HARQ indicator
            convert_to_bits(0,&dci,3);//Only considering FDD
            //New data indicator
            convert_to_bits(alloc->ndi,&dci,1);
            //Redundancy version
            convert_to_bits(alloc->rv_idx,&dci,2);
            //TPC
            convert_to_bits(alloc->tpc,&dci,2);

            alloc->tbs=TBS[alloc->mcs][alloc->N_prb-1];
        }

    size =dci-out_bits;
    if(size == 12 ||
    size == 14 ||
    size == 16 ||
    size == 20 ||
    size == 24 ||
    size == 26 ||
    size == 32 ||
    size == 40 ||
    size == 44 ||
    size == 56)
    {size++;
        convert_to_bits(0,&dci,1);
    }
    *N_out_bits=size;
}
void dci_0_pack(PHY_ALLOCATION_STRUCT *alloc,int N_rb_ul,int N_ant,int *out_bits,int *N_out_bits)
{
    int  RIV;
    int  RIV_length;
    int  size;
    int  *dci = out_bits;

    // Carrier indicator

        convert_to_bits(0, &dci, 3);

    // Format 0/1A flag is set to format 0
    convert_to_bits(DCI_0_1A_FLAG_0, &dci, 1);

    // Frequency hopping flag
    convert_to_bits(0, &dci, 1);

    // RBA
    // FIXME: Only supporting non-hopping single-cluster
    RIV_length = ceilf(logf(N_rb_ul*(N_rb_ul+1)/2)/logf(2));
    if((alloc->N_prb-1) <= (N_rb_ul/2))
    {
        RIV = N_rb_ul*(alloc->N_prb-1) + alloc->prb[0][0];
    }else{
        RIV = N_rb_ul*(N_rb_ul - alloc->N_prb + 1) + (N_rb_ul - 1 - alloc->prb[0][0]);
    }
    convert_to_bits(RIV, &dci, RIV_length);

    // Modulation and coding scheme and redundancy version
    convert_to_bits(alloc->mcs, &dci, 5);

    // New data indicator
    convert_to_bits(alloc->ndi, &dci, 1);

    // TPC command
    convert_to_bits(alloc->tpc, &dci, 2);

    // Cyclic shift
    convert_to_bits(0, &dci, 3);

    // CSI request
    convert_to_bits(0, &dci, 1);

    // Pad bit
    convert_to_bits(0, &dci, 1);

    // Pad if needed
    size = dci - out_bits;
    if(size == 12 ||
       size == 14 ||
       size == 16 ||
       size == 20 ||
       size == 24 ||
       size == 26 ||
       size == 32 ||
       size == 40 ||
       size == 44 ||
       size == 56)
    {
        size++;
        convert_to_bits(0, &dci, 1);
    }
    *N_out_bits = size;
}

void calculate_cce(PHY_PCFICH *pcfich,PHY_PHICH *phich,LIBLTE_PHY *phy_struct,int N_ant,int pdcch_symbols,int *N_cce,int *N_reg_pdcch)
    {
        int N_phich_reg;
        int N_pcfich_reg=pcfich->N_reg;
        int N_reg_rb=3;
        N_phich_reg=phy_struct->N_group_phich*3;
        *N_reg_pdcch=pdcch_symbols*(phy_struct->N_rb_dl*N_reg_rb)-phy_struct->N_rb_dl-N_phich_reg-N_pcfich_reg;

        if(N_ant==4)
            *N_reg_pdcch-=phy_struct->N_rb_dl;
            *N_cce = (*N_reg_pdcch/N_REG_CCE);
    }

void dci_encode(LIBLTE_PHY *phy_struct,int *in_bits,int rnti,int ue_ant,int N_in_bits,int *out_bits,int N_out_bits)
{  int *a_bits;
    int i,N_d_bits;
    int p_bits[16];
    int x_rnti[16];
    int x_as_bits[16];
    int g[3] = {0133, 0171, 0165};
    memset(x_as_bits,0,sizeof(int)*16);
    a_bits=in_bits;
    calc_crc(&a_bits[0],N_in_bits,CRC16,&p_bits[0],16);
    for(i=0;i<16;i++)
        x_rnti[i]=(rnti>>(16-i))&0x1;
    if(ue_ant)
        x_as_bits[15]=1;

    //Masking the sequence
    for(i=0;i<16;i++)
        p_bits[i]^=x_as_bits[i]^x_rnti[i];
    for(i=0;i<N_in_bits;i++)
        phy_struct->dci_c[i]=in_bits[i];
    for(i=0;i<16;i++)
        phy_struct->dci_c[N_in_bits+i]=p_bits[i];


conv_encode(phy_struct,phy_struct->dci_c,N_in_bits+16,7,3,g,1,phy_struct->dci_tx_d_bits,&N_d_bits);

    rate_match_conv(phy_struct,phy_struct->dci_tx_d_bits,N_d_bits,N_out_bits,&out_bits[0]);

}

void pdcch_pre_calc(LIBLTE_PHY *phy_struct,int N_ant)
{
    int N_reg_phich,N_reg_pdcch;
    int N_symbs;
    int i,p;
    int N_reg_pcfich=4;

    N_reg_phich=phy_struct->N_group_phich*3;
    for(N_symbs=1;N_symbs<3;N_symbs++)
    {
        N_reg_pdcch=N_symbs*(3*phy_struct->N_rb_dl)-(phy_struct->N_rb_dl)-(N_reg_phich)-(N_reg_pcfich);
        for(i=0;i<N_reg_pdcch;i++)
        {
            phy_struct->pdcch_reg_vec[i]=i;
        }
        int R_c;
        int C_c=32;
        int i,j,k,x,N_dummy;
        int i_dx;
        int w_idx=0;
        int K_pi,K_w;
        R_c=0;

        while(N_reg_pdcch>R_c*C_c)
        {
            R_c++;
        }
        for(x=0;x<3;x++)
        {
        if(R_c*C_c>N_reg_pdcch)
            N_dummy=R_c*C_c-N_reg_pdcch;
        else
            N_dummy=0;
        for(i=0;i<N_dummy;i++)
            phy_struct->rmc_tmp[i]=0;
        i_dx=0;
        for(i=N_dummy;i<C_c*R_c;i++)
        {
            phy_struct->rmc_tmp[i]=phy_struct->pdcch_reg_vec[i];
            i_dx++;
        }



            i_dx=0;
        for(i=0;i<R_c;i++)
        { for(j=0;j<C_c;j++)
          {
            phy_struct->rmc_sb_mat[i][j]=phy_struct->rmc_tmp[i_dx++];


           }

        }

        //Permutation
        for(i=0;i<R_c;i++)
        {
            for(j=0;j<C_c;j++)
            {
                phy_struct->rmc_sb_perm_mat[i][j]=phy_struct->rmc_sb_mat[i][ICC_PERM[j]];

            }

        }
        for(j=0;j<C_c;j++)
        {
            for(i=0;i<R_c;i++)
            {
            phy_struct->rmc_w[w_idx++]=phy_struct->rmc_sb_perm_mat[i][j];

            }


            }


    }
     //Creating the circular buffer
        K_pi = R_c*C_c;
        k=0;
        j=0;
        while(k<N_reg_pdcch)
        {
            if(phy_struct->rmc_w[j%K_pi]!=0)
            {
                phy_struct->pdcch_reg_perm_vec[k++]=phy_struct->rmc_w[j%K_pi];


            }
            j++;
        }

        for(p=0;p<N_ant;p++)
        {
            for(i=0;i<N_reg_pdcch;i++)
                phy_struct->pdcch_permute_map[N_reg_pdcch][i]=phy_struct->pdcch_reg_perm_vec[i];
        }

    }
}

void pdcch_channel_encode(PHY_PCFICH *pcfich,PHY_PHICH *phich,LIBLTE_PHY *phy_struct,PDCCH_STRUCT *pdcch,FRAME *subframe,int N_ant,int N_id_cell)
{   int c_init,a_idx,dci_size,N_bits,c_idx,M_symb,M_layer_symb,M_ap_symb,p,i_dx;
    int pdcch_symbols,i,j,Y_k,actual_idx,a_dx,k;
    int shift_idx;
    int N_cce;
    int N_reg_pdcch;
    int k_prime,l_prime,m_prime;
    bool valid_reg;
    pdcch_symbols=pcfich->cfi;

    calculate_cce(pcfich,phich,phy_struct,N_ant,pdcch_symbols,&N_cce,&N_reg_pdcch);
    for(i=0;i<N_cce;i++)
    {
      for(j=0;j<N_REG_CCE*4;j++)
       {
        phy_struct->pdcch_cce_re[i][j]=0;
        phy_struct->pdcch_cce_im[i][j]=0;


       }
        phy_struct->pdcch_used[i]=false;

    }

    c_init=(subframe->num<<9)+N_id_cell;
    scramb_generate(c_init,PDCCH_BITS_MAX*2,phy_struct->pdcch_c); //Generates a scrambling sequence of length equal to max of the search space i.e 16CCE's

    for(a_idx=0;a_idx<pdcch->N_alloc;a_idx++)
    {
        if(pdcch->alloc[a_idx].chan_type==DLSCH)
            dci_1a_pack(&pdcch->alloc[a_idx],phy_struct->N_rb_dl,N_ant,phy_struct->pdcch_dci,&dci_size);

        else
            dci_0_pack(&pdcch->alloc[a_idx],phy_struct->N_rb_ul,N_ant,phy_struct->pdcch_dci,&dci_size);


        N_bits=288;

        dci_encode(phy_struct,phy_struct->pdcch_dci,pdcch->alloc[a_idx].rnti,0,dci_size,phy_struct->pdcch_encode_bits,N_bits);

        //check for user specific and common search space mapping and RE mapping

        //Common search space and user specific search space

        if(MAC_PRNTI==pdcch->alloc[a_idx].rnti||MAC_SIRNTI==pdcch->alloc[a_idx].rnti||(MAC_SI_RNTI_START<=pdcch->alloc[a_idx].rnti&&MAC_SI_RNTI_END>=pdcch->alloc[a_idx].rnti))    //Common search space
        {
            for(c_idx=0;c_idx<4;c_idx++)
            {

                if(!phy_struct->pdcch_used[4*c_idx+0]&&!phy_struct->pdcch_used[4*c_idx+1]&&!phy_struct->pdcch_used[4*c_idx+2]&&!phy_struct->pdcch_used[4*c_idx+3])
                {
                    for(i=0;i<N_bits;i++)
                        phy_struct->pdcch_scramb_bits[i]=phy_struct->pdcch_encode_bits[i]^phy_struct->pdcch_c[4*c_idx*N_REG_CCE*4*2+i];

                    modulation_mapper(phy_struct->pdcch_c,N_bits,MODULATION_TYPE_QPSK,phy_struct->pdcch_d_re,phy_struct->pdcch_d_im,&M_symb);

               layer_mapper(phy_struct->pdcch_d_re,phy_struct->pdcch_d_im,M_symb,N_ant,1,transmit,phy_struct->pdcch_x_re,phy_struct->pdcch_x_im,&M_layer_symb);

                pre_coder_dl(phy_struct->pdcch_x_re,phy_struct->pdcch_x_im,M_layer_symb,N_ant,transmit,phy_struct->pdcch_y_re[0],phy_struct->pdcch_y_im[0],576,&M_ap_symb);


                    //Mapping to  the CCE ,RE locations

                    for(p=0;p<N_ant;p++)
                    {  i_dx=0;
                        for(i=0;i<4;i++)
                          {
                              for(j=0;j<N_REG_CCE*4;j++)
                              {
                                  phy_struct->pdcch_re[p][4*c_idx+i][j]=phy_struct->pdcch_y_re[p][i_dx];
                                  phy_struct->pdcch_im[p][4*c_idx+i][j]=phy_struct->pdcch_y_im[p][i_dx];

                                  i_dx++;
                              }

                              phy_struct->pdcch_used[4*c_idx+i]=true;

                            }

                }
                    break;
            }




        }


    } else  //User search space
    {    Y_k=pdcch->alloc[a_idx].rnti;
         for(i=0;i<subframe->num;i++)
         {
             Y_k=(39827*Y_k)%65537;
         }

         for(a_dx=0;a_dx<2;i++)
         {
             actual_idx=4*((Y_k+i)%(N_cce/4));

             if(!phy_struct->pdcch_used[actual_idx+0]&&!phy_struct->pdcch_used[actual_idx+1]&&!phy_struct->pdcch_used[actual_idx+2]&&!phy_struct->pdcch_used[actual_idx+3])
             {

                 for(i=0;i<N_bits;i++)
                 phy_struct->pdcch_scramb_bits[i]=phy_struct->pdcch_encode_bits[i]^phy_struct->pdcch_c[4*a_dx*N_REG_CCE*4*2+i];

            modulation_mapper(phy_struct->pdcch_c,N_bits,MODULATION_TYPE_QPSK,phy_struct->pdcch_d_re,phy_struct->pdcch_d_im,&M_symb);

            layer_mapper(phy_struct->pdcch_d_re,phy_struct->pdcch_d_im,M_symb,N_ant,1,transmit,phy_struct->pdcch_x_re,phy_struct->pdcch_x_im,&M_layer_symb);

            pre_coder_dl(phy_struct->pdcch_x_re,phy_struct->pdcch_x_im,M_layer_symb,N_ant,transmit,phy_struct->pdcch_y_re[0],phy_struct->pdcch_y_im[0],576,&M_ap_symb);


                //Mapping to  the CCE ,RE locations

                for(p=0;p<N_ant;p++)
                  {  i_dx=0;
                     for(i=0;i<4;i++)
                        {
                            for(j=0;j<N_REG_CCE*4;j++)
                                {
                                    phy_struct->pdcch_re[p][actual_idx+i][j]=phy_struct->pdcch_y_re[p][i_dx];
                                    phy_struct->pdcch_im[p][actual_idx+i][j]=phy_struct->pdcch_y_im[p][i_dx];

                                    i_dx++;
                                }

                            phy_struct->pdcch_used[actual_idx+i]=true;

                        }

                    }
                 break;
             }

         }

    }

        //RE mapping

        //Group into quadrapulets
        for(p=0;p<N_ant;p++)
        {
            for(i=0;i<N_cce;i++)
            {
                for(j=0;j<N_reg_pdcch;j++)
                {
                    for(k=0;k<4;k++)
                {
                    phy_struct->pdcch_reg_re[p][i*N_REG_CCE+j][k]=phy_struct->pdcch_re[p][i][4*j+k];
                }
                }


            }
        }

        //Permutation

        for(p=0;p<N_ant;p++)
        {
            for(i=0;i<N_reg_pdcch;i++)
            {
                for(j=0;j<4;j++)
                {
                    phy_struct->pdcch_perm_re[p][i][j]=phy_struct->pdcch_reg_re[p][phy_struct->pdcch_permute_map[N_reg_pdcch][i]][j];
                    phy_struct->pdcch_perm_im[p][i][j]=phy_struct->pdcch_reg_im[p][phy_struct->pdcch_permute_map[N_reg_pdcch][i]][j];

                }
            }
        }


}
    //Cyclically shifting

    for(p=0;p<N_ant;p++)
    {
        for(i=0;i<N_reg_pdcch;i++)
        {    shift_idx=(i+N_id_cell)%N_reg_pdcch;
            for(j=0;j<4;j++)
            {
                phy_struct->pdcch_shift_re[p][i][j]=phy_struct->pdcch_perm_re[p][shift_idx][j];
                phy_struct->pdcch_shift_im[p][i][j]=phy_struct->pdcch_perm_im[p][shift_idx][j];

            }

        }
    }

 //REG mapping
    m_prime=0;
    k_prime=0;
    while(k_prime<phy_struct->N_rb_dl*phy_struct->N_sc_rb)
    {
        l_prime=0;
        valid_reg=true;
        while(l_prime<pdcch->N_symbs)
        {
         if(l_prime==0)
          {
            for(i=0;i<pcfich->N_reg;i++)
            {
            if(k_prime==pcfich->k[i])
                valid_reg=false;
            }

            for(i=0;i<phich->N_reg;i++)
            {
                if(k_prime==phich->k[i])
                    valid_reg=false;
            }

            if(valid_reg==true&&(k_prime%6==0))
            {

                if(m_prime<N_reg_pdcch)
                {
                    for(p=0;p<N_ant;p++)
                    {  i_dx=0;
                        for(i=0;i<6;i++)
                        {   //Avoid CRS
                            if((N_id_cell)%3!=i%3)
                            {
                                subframe->tx_symb_re[p][l_prime][k_prime+i]=phy_struct->pdcch_shift_re[p][m_prime][i_dx];
                                i_dx++;

                            }

                        }
                    }
                    m_prime++;
                }

            }


        }

        else if(l_prime==1&&N_ant==4)
          {
            if(k_prime%6==0)
              {  if(m_prime<N_reg_pdcch)
               {
                  for(p=0;p<N_ant;p++)
                  {   i_dx=0;
                      for(i=0;i<6;i++)
                      {
                          if((N_id_cell)%3!=i%3)
                          {
                            subframe->tx_symb_re[p][l_prime][k_prime+i]=phy_struct->pdcch_shift_re[p][m_prime][i_dx];
                              i_dx++;

                          }
                      }
                  }

                }
                  m_prime++;
              }
          }
        else
        {
            if(k_prime%4==0)
            {  if(m_prime<N_reg_pdcch)
             {
                for(p=0;p<N_ant;p++)
                {   i_dx=0;
                    for(i=0;i<6;i++)
                    {
                        if((N_id_cell)%3!=i%3)
                        {
                          subframe->tx_symb_re[p][l_prime][k_prime+i]=phy_struct->pdcch_shift_re[p][m_prime][i_dx];
                            i_dx++;

                        }
                    }
                }

              }
                m_prime++;
            }
        }
            l_prime++;
    }

        k_prime++;






    }

}



void generate_dl_subframe(LIBLTE_PHY *phy_struct,FRAME *subframe,int ant,float *i_samps,float*q_samps)
{
    int i;
    int N_samps;
    int idx=0;
    for(i=0;i<14;i++)
    {
        idx+=N_samps;
        symbols_to_samples(phy_struct,&subframe->tx_symb_re[ant][i][0],&subframe->tx_symb_re[ant][i][0],i,&i_samps[idx],&q_samps[idx],&N_samps);
    }

}
void wrap_phase(float *phase1,float phase2)
{
    while((*phase1-phase2)>M_PI)

        *phase1-=2*M_PI;

    while((*phase1-phase2)<-M_PI)
        *phase1+=2*M_PI;
}

void get_dl_subframe_and_ce(LIBLTE_PHY *phy_struct,float *i_samps,float *q_samps,int frame_start_idx,int subfr_num,int N_id_cell,int N_ant,FRAME *subframe)
{
    int i,j,z;
    int v[5];
    float *symb_re,*symb_im,*rs_re,*rs_im;
    int N_symb;
    float ce_mag,ce_ang;
    int sym[5];
    int k,m_prime;
    float tmp_re,tmp_im;
    float frac_mag,frac_phase;
    subframe->num=subfr_num;
    int v_shift=N_id_cell%6;
    int subfr_start_index=frame_start_idx+subframe->num*phy_struct->N_samps_per_subfr;

    for(i=0;i<16;i++)
    {
        samples_to_symbols(phy_struct, i_samps,q_samps, i%7,subfr_start_index+(i/7)*phy_struct->N_samps_per_slot, &subframe->rx_symb_re[i][0], &subframe->rx_symb_im[i][0]);

    }
    generate_crs((subframe->num*2+0)%20,0,0,N_id_cell,phy_struct->dl_ce_crs_re[0],phy_struct->dl_ce_crs_im[0]);
    generate_crs((subframe->num*2+0)%20,1,0,N_id_cell,phy_struct->dl_ce_crs_re[1],phy_struct->dl_ce_crs_im[1]);
    generate_crs((subframe->num*2+0)%20,4,0,N_id_cell,phy_struct->dl_ce_crs_re[4],phy_struct->dl_ce_crs_im[4]);
    generate_crs((subframe->num*2+1)%20,0,0,N_id_cell,phy_struct->dl_ce_crs_re[7],phy_struct->dl_ce_crs_im[7]);
    generate_crs((subframe->num*2+1)%20,1,0,N_id_cell,phy_struct->dl_ce_crs_re[8],phy_struct->dl_ce_crs_im[8]);
    generate_crs((subframe->num*2+1)%20,4,0,N_id_cell,phy_struct->dl_ce_crs_re[11],phy_struct->dl_ce_crs_im[11]);
    generate_crs((subframe->num*2+2)%20,0,0,N_id_cell,phy_struct->dl_ce_crs_re[14],phy_struct->dl_ce_crs_im[14]);
    generate_crs((subframe->num*2+2)%20,1,0,N_id_cell,phy_struct->dl_ce_crs_re[15],phy_struct->dl_ce_crs_im[15]);

    if(N_ant==0)
    {
        v[0]=0;
        v[1]=3;
        v[2]=0;
        v[3]=3;
        v[4]=0;
        sym[0]=0;
        sym[1]=4;
        sym[2]=7;
        sym[3]=11;
        sym[4]=14;
        N_symb=5;

    }
    else if (N_ant==1)
    {
        v[0]=3;
        v[1]=0;
        v[2]=3;
        v[3]=0;
        v[4]=3;
        sym[0]=0;
        sym[1]=4;
        sym[2]=7;
        sym[3]=11;
        sym[4]=14;
        N_symb=5;
    }
    else if (N_ant==2)
    {
        v[0]=0;
        v[1]=3;
        v[2]=0;
        sym[0]=1;
        sym[1]=8;
        sym[2]=15;
        N_symb=3;
    }
    else if (N_ant==3)
    {
        v[0]   = 3;
        v[1]   = 6;
        v[2]   = 3;
        sym[0] = 1;
        sym[1] = 8;
        sym[2] = 15;
        N_symb  = 3;

    }
    for(i=0;i<N_symb;i++)
    {
        symb_re=&subframe->rx_symb_re[sym[i]][0];
        symb_im=&subframe->rx_symb_re[sym[i]][0];
        rs_re=&phy_struct->dl_ce_crs_re[sym[i]][0];
        rs_im=&phy_struct->dl_ce_crs_im[sym[i]][0];

        for(j=0;j<2*phy_struct->N_rb_dl;i++)
        {
            k=6*j+(v[i]+v_shift)%6;
            m_prime = j + N_rb_max - phy_struct->N_rb_dl;

            tmp_re=symb_re[k]*rs_re[m_prime]+symb_im[k]*rs_im[m_prime];
            tmp_im=symb_re[k]*rs_im[m_prime]+symb_im[k]*rs_re[m_prime];
            phy_struct->dl_ce_mag[i][k]=sqrt(tmp_re*tmp_re+tmp_im*tmp_im);
            phy_struct->dl_ce_ang[i][k]=atan2f(tmp_im,tmp_re);

            if(j>0)
            {
                wrap_phase(&phy_struct->dl_ce_ang[i][k],phy_struct->dl_ce_ang[i][k-6]);
                //Interpolate between the CRS
                frac_mag=(phy_struct->dl_ce_mag[i][k]-phy_struct->dl_ce_mag[i][k-6])/6;
                frac_phase=(phy_struct->dl_ce_ang[i][k]-phy_struct->dl_ce_ang[i][k-6])/6;

                for(z=1;z<6;z++)
                {
                    phy_struct->dl_ce_mag[i][k-z]=phy_struct->dl_ce_mag[i][k-(z-1)]-frac_mag;
                    phy_struct->dl_ce_ang[i][k-z]=phy_struct->dl_ce_ang[i][k-(z-1)]-frac_mag;
                }


            }

            if(j==1)
            {
                for(z=1;z<(v[i]+v_shift)%6+1;z++)
                {
                    phy_struct->dl_ce_mag[i][k-6-z]=phy_struct->dl_ce_mag[i][k-6-(z-1)]-frac_mag;
                    phy_struct->dl_ce_ang[i][k-6-z]=phy_struct->dl_ce_ang[i][k-6-(z-1)]-frac_mag;
                }


            }
        }

        for(z=1;z<5-(v[i]+v_shift)%6+1;z++)
        {
            phy_struct->dl_ce_mag[i][k+z]=phy_struct->dl_ce_mag[i][k+(z-1)]-frac_mag;
            phy_struct->dl_ce_ang[i][k+z]=phy_struct->dl_ce_ang[i][k+(z-1)]-frac_mag;
        }


    }
    //Linearly interpolate betweeen the symbols

    if(N_symb==3) //Directly get the channel eestimates for symbols 1 and 8
    {
        for(j=0;j<phy_struct->N_rb_dl*phy_struct->N_sc_rb;j++)
        {
            subframe->rx_ce_re[N_ant][1][j]=phy_struct->dl_ce_mag[0][j]*cosf(phy_struct->dl_ce_ang[0][j]);
            subframe->rx_ce_im[N_ant][1][j]=phy_struct->dl_ce_mag[0][j]*sinf(phy_struct->dl_ce_ang[0][j]);
            subframe->rx_ce_re[N_ant][8][j]=phy_struct->dl_ce_mag[1][j]*cosf(phy_struct->dl_ce_ang[1][j]);
            subframe->rx_ce_im[N_ant][8][j]=phy_struct->dl_ce_mag[1][j]*sinf(phy_struct->dl_ce_ang[1][j]);
            //Interpolate for 2,3,4,5,6,7
            frac_mag=(phy_struct->dl_ce_mag[1][j]-phy_struct->dl_ce_mag[0][j])/7;
            wrap_phase(&phy_struct->dl_ce_mag[1][j], phy_struct->dl_ce_mag[0][j]);
            frac_phase=(phy_struct->dl_ce_mag[1][j]-phy_struct->dl_ce_mag[0][j]);
            wrap_phase(&frac_phase,0);
            frac_phase/=7;
            ce_mag=phy_struct->dl_ce_mag[1][j];
            ce_ang=phy_struct->dl_ce_ang[1][j];
            for(z=7;z>1;z++)
            {   ce_mag-=frac_mag;
                ce_ang-=frac_phase;
                subframe->rx_ce_re[N_ant][z][j]=ce_mag*cosf(ce_ang);
                subframe->rx_ce_im[N_ant][z][j]=ce_mag*sinf(ce_ang);

            }
            //Interpolating for symbol 0
            //FIX ME:should ideally have used the previous slot for finding the estimate
            ce_mag=phy_struct->dl_ce_mag[0][j]-frac_mag;
            ce_ang=phy_struct->dl_ce_ang[0][j]-frac_phase;
            subframe->rx_ce_re[N_ant][0][j]=ce_mag*cosf(ce_ang);
            subframe->rx_ce_im[N_ant][0][j]=ce_mag*sinf(ce_ang);
            //Interpolating for 9,10,11,12,13

            frac_mag=(phy_struct->dl_ce_mag[2][j]-phy_struct->dl_ce_mag[1][j])/7;
            wrap_phase(&phy_struct->dl_ce_mag[2][j], phy_struct->dl_ce_mag[1][j]);
            frac_phase=(phy_struct->dl_ce_mag[2][j]-phy_struct->dl_ce_mag[1][j]);
            wrap_phase(&frac_phase,0);
            frac_phase/=7;
            ce_mag=phy_struct->dl_ce_mag[2][j]-frac_mag;
            ce_ang=phy_struct->dl_ce_ang[2][j]-frac_phase;
            for(z=13;z>8;z++)
            {   ce_mag-=frac_mag;
                ce_ang-=frac_phase;
                subframe->rx_ce_re[N_ant][z][j]=ce_mag*cosf(ce_ang);
                subframe->rx_ce_im[N_ant][z][j]=ce_mag*sinf(ce_ang);

            }

        }


    }
    else  //Directly assign for 0,4,7,11
    {
        for(j=0;j<phy_struct->N_rb_dl*phy_struct->N_sc_rb;j++)
        {
            subframe->rx_ce_re[N_ant][0][j]=phy_struct->dl_ce_mag[0][j]*cosf(phy_struct->dl_ce_ang[0][j]);
            subframe->rx_ce_im[N_ant][0][j]=phy_struct->dl_ce_mag[0][j]*sinf(phy_struct->dl_ce_ang[0][j]);
            subframe->rx_ce_re[N_ant][4][j]=phy_struct->dl_ce_mag[1][j]*cosf(phy_struct->dl_ce_ang[1][j]);
            subframe->rx_ce_im[N_ant][4][j]=phy_struct->dl_ce_mag[1][j]*sinf(phy_struct->dl_ce_ang[1][j]);
            subframe->rx_ce_re[N_ant][7][j]=phy_struct->dl_ce_mag[2][j]*cosf(phy_struct->dl_ce_ang[2][j]);
            subframe->rx_ce_im[N_ant][7][j]=phy_struct->dl_ce_mag[2][j]*sinf(phy_struct->dl_ce_ang[2][j]);
            subframe->rx_ce_re[N_ant][11][j]=phy_struct->dl_ce_mag[3][j]*cosf(phy_struct->dl_ce_ang[3][j]);
            subframe->rx_ce_im[N_ant][11][j]=phy_struct->dl_ce_mag[3][j]*sinf(phy_struct->dl_ce_ang[3][j]);

           //Interpolate for 1,2,3,5,6,8,9,10
           frac_mag=(phy_struct->dl_ce_mag[1][j]-phy_struct->dl_ce_mag[0][j])/4;
           wrap_phase(&phy_struct->dl_ce_mag[1][j], phy_struct->dl_ce_mag[0][j]);
           frac_phase=(phy_struct->dl_ce_mag[1][j]-phy_struct->dl_ce_mag[0][j]);
            wrap_phase(&frac_phase,0);
            frac_phase/=4;
            ce_mag=phy_struct->dl_ce_mag[1][j];
            ce_ang=phy_struct->dl_ce_ang[1][j];
            for(z=3;z>0;z++)
            {   ce_mag-=frac_mag;
                ce_ang-=frac_phase;
                subframe->rx_ce_re[N_ant][z][j]=ce_mag*cosf(ce_ang);
                subframe->rx_ce_im[N_ant][z][j]=ce_mag*sinf(ce_ang);

            }
            frac_mag=(phy_struct->dl_ce_mag[2][j]-phy_struct->dl_ce_mag[1][j])/3;
            wrap_phase(&phy_struct->dl_ce_mag[2][j], phy_struct->dl_ce_mag[1][j]);
            frac_phase=(phy_struct->dl_ce_mag[2][j]-phy_struct->dl_ce_mag[1][j]);
             wrap_phase(&frac_phase,0);
             frac_phase/=3;
             ce_mag=phy_struct->dl_ce_mag[2][j];
             ce_ang=phy_struct->dl_ce_ang[2][j];
             for(z=6;z>4;z++)
             {   ce_mag-=frac_mag;
                 ce_ang-=frac_phase;
                 subframe->rx_ce_re[N_ant][z][j]=ce_mag*cosf(ce_ang);
                 subframe->rx_ce_im[N_ant][z][j]=ce_mag*sinf(ce_ang);

             }
            frac_mag=(phy_struct->dl_ce_mag[3][j]-phy_struct->dl_ce_mag[2][j])/4;
            wrap_phase(&phy_struct->dl_ce_mag[3][j], phy_struct->dl_ce_mag[2][j]);
            frac_phase=(phy_struct->dl_ce_mag[3][j]-phy_struct->dl_ce_mag[2][j]);
             wrap_phase(&frac_phase,0);
             frac_phase/=4;
             ce_mag=phy_struct->dl_ce_mag[3][j];
             ce_ang=phy_struct->dl_ce_ang[3][j];
             for(z=10;z>7;z++)
             {   ce_mag-=frac_mag;
                 ce_ang-=frac_phase;
                 subframe->rx_ce_re[N_ant][z][j]=ce_mag*cosf(ce_ang);
                 subframe->rx_ce_im[N_ant][z][j]=ce_mag*sinf(ce_ang);

             }
            //Interpolate for 12,13
            frac_mag=(phy_struct->dl_ce_mag[4][j]-phy_struct->dl_ce_mag[3][j])/4;
            wrap_phase(&phy_struct->dl_ce_mag[4][j], phy_struct->dl_ce_mag[3][j]);
            frac_phase=(phy_struct->dl_ce_mag[4][j]-phy_struct->dl_ce_mag[3][j]);
             wrap_phase(&frac_phase,0);
             frac_phase/=4;
             ce_mag=phy_struct->dl_ce_mag[4][j]-frac_mag;
             ce_ang=phy_struct->dl_ce_ang[4][j]-frac_phase;
             for(z=13;z>11;z++)
             {   ce_mag-=frac_mag;
                 ce_ang-=frac_phase;
                 subframe->rx_ce_re[N_ant][z][j]=ce_mag*cosf(ce_ang);
                 subframe->rx_ce_im[N_ant][z][j]=ce_mag*sinf(ce_ang);

             }



        }
    }





}
void pcfich_channel_decode(LIBLTE_PHY *phy_struct,int N_id_cell,int N_ant,PHY_PCFICH *pcfich,FRAME *subframe,int *N_bits)
{
    int k_bar,i,p,idx,j;
    int c_init;
    pcfich->N_reg=4;
    k_bar=(phy_struct->N_sc_rb/2)*(N_id_cell%(2*phy_struct->N_rb_dl));

       for(i=0;i<pcfich->N_reg;i++)
       {
pcfich->k[i]=(k_bar+(i*phy_struct->N_rb_dl/2)*(phy_struct->N_sc_rb/2))%((phy_struct->N_sc_rb)*(phy_struct->N_rb_dl));

           pcfich->n[i]=(pcfich->k[i]/6)-0.5;
           idx=0;

               for(j=0;j<6;j++)
               {
                   if((N_id_cell%3)!=(j%3))
                   {
                       phy_struct->pdcch_y_re_est[4*i+idx]=subframe->rx_symb_re[0][pcfich->k[i]+j];
                       phy_struct->pdcch_y_im_est[4*i+idx]=subframe->rx_symb_im[0][pcfich->k[i]+j];
                       for(p=0;p<N_ant;p++)
                       {
                           phy_struct->pdcch_c_re_est[p][4*i+idx]=subframe->rx_ce_re[p][0][pcfich->k[i]+j];
                           phy_struct->pdcch_c_im_est[p][4*i+idx]=subframe->rx_ce_im[p][0][pcfich->k[i]+j];
                       }
                       idx++;
                   }
               }

       }

    c_init=((subframe->num)+1)*((2*N_id_cell+1)<<9) +N_id_cell;
    scramb_generate(c_init,32,phy_struct->pdcch_c);





}

void viterbi_decode(LIBLTE_PHY *phy_struct,
                    int             *d_bits,
                    int             N_d_bits,
                    int             constraint_len,
                    int             rate,
                    int            *g,
                    int             *c_bits,
                    int            *N_c_bits)
{
    float  init_min;
    int len2;
    float  tmp1;
    float  tmp2;
    int  i;
    int j;
    int k;
    int o;
    int N_states = 1<<(constraint_len-1);
    int idx;
    int tb_state_len;
    int  in_path;
    int  prev_state;
    int  in_bit;
    int  prev_state_0;
    int  prev_state_1;
    int  s_reg[constraint_len];
    int  g_array[3][constraint_len];

    // Convert g to binary
    for(i=0; i<rate; i++)
    {
        for(j=0; j<constraint_len; j++)
        {
            g_array[i][j] = (g[i] >> (constraint_len-j-1)) & 1;
        }
    }

    // Precalculate state transition outputs
    for(i=0; i<N_states; i++)
    {
        // Determine the input path
        if(i < (N_states/2))
        {
            in_path = 0;
        }else{
            in_path = 1;
        }

        // Determine the outputs based on the previous state and input path
        for(j=0; j<2; j++)
        {
            prev_state = ((i << 1) + j) % N_states;
            for(k=0; k<constraint_len; k++)
            {
                s_reg[k] = (prev_state >> (constraint_len-k-1)) & 1;
            }
            s_reg[0] = in_path;
            for(k=0; k<rate; k++)
            {
                phy_struct->vd_st_output[i][j][k] = 0;
                for(o=0; o<constraint_len; o++)
                {
                    phy_struct->vd_st_output[i][j][k] += s_reg[o]*g_array[k][o];
                }
                phy_struct->vd_st_output[i][j][k] %= 2;

            }
        }
    }

    // Calculate branch and path metrics
    for(i=0; i<N_states; i++)
    {
        for(j=0; j<(N_d_bits/rate); j++)
        {
            phy_struct->vd_path_metric[i][j] = 0;
        }
    }
    for(i=0; i<(N_d_bits/rate); i++)
    {
        for(j=0; j<N_states; j++)
        {
            phy_struct->vd_br_metric[j][0] = 0;
            phy_struct->vd_br_metric[j][1] = 0;
            phy_struct->vd_p_metric[j][0]  = 0;
            phy_struct->vd_p_metric[j][1]  = 0;
            phy_struct->vd_w_metric[j][0]  = 0;
            phy_struct->vd_w_metric[j][1]  = 0;

            // Calculate the accumulated branch metrics for each state
            for(k=0; k<2; k++)
            {
                prev_state                    = ((j<<1)+k) % N_states;
                phy_struct->vd_p_metric[j][k] = phy_struct->vd_path_metric[prev_state][i];
                for(o=0; o<rate; o++)
                {
                    if(d_bits[i*rate + o] > 0)
                    {
                        in_bit = 1;
                    }else{
                        in_bit = 0;
                    }
                    phy_struct->vd_br_metric[j][k] += (phy_struct->vd_st_output[j][k][o]+in_bit)%2;

                    phy_struct->vd_w_metric[j][k]  += fabs(d_bits[i*rate + o]);
                }

            }

            // Keep the smallest branch metric as the path metric, weight the branch metric
            tmp1 = phy_struct->vd_br_metric[j][0] + phy_struct->vd_p_metric[j][0];
            tmp2 = phy_struct->vd_br_metric[j][1] + phy_struct->vd_p_metric[j][1];
            if(tmp1 > tmp2)
            {
                phy_struct->vd_path_metric[j][i+1] = phy_struct->vd_p_metric[j][1] + phy_struct->vd_w_metric[j][1]*phy_struct->vd_br_metric[j][1];
            }else{
                phy_struct->vd_path_metric[j][i+1] = phy_struct->vd_p_metric[j][0] + phy_struct->vd_w_metric[j][0]*phy_struct->vd_br_metric[j][0];
            }

        }
    }

    // Find the minimum metric for the last iteration
    init_min                     = 1000000;
    idx                          = 0;
    phy_struct->vd_tb_state[idx] = 1000000;
    for(i=0; i<N_states; i++)
    {
        if(phy_struct->vd_path_metric[i][(N_d_bits/rate)] < init_min)
        {
            init_min                       = phy_struct->vd_path_metric[i][(N_d_bits/rate)];
            phy_struct->vd_tb_state[idx++] = i;
        }
    }
    len2=idx-1;
    // Traceback to find the minimum path metrics at each iteration
    for(i=(N_d_bits/rate)-1; i>=0; i--)
    {
        prev_state_0 = ((((int)phy_struct->vd_tb_state[idx-1])<<1) + 0) % N_states;
        prev_state_1 = ((((int)phy_struct->vd_tb_state[idx-1])<<1) + 1) % N_states;

        // Keep the smallest state
        if(phy_struct->vd_path_metric[prev_state_0][i] > phy_struct->vd_path_metric[prev_state_1][i])
        {
            phy_struct->vd_tb_state[idx++] = prev_state_1;
        }else{
            phy_struct->vd_tb_state[idx++] = prev_state_0;
        }
    }
    tb_state_len = idx;
    // Read through the traceback to determine the input bits
    idx = 0;
    for(i=tb_state_len-2; i>=len2; i--)
    {
        // If transition has resulted in a lower valued state,
        // the output is 0 and vice-versa
        if(phy_struct->vd_tb_state[i] < phy_struct->vd_tb_state[i+1])
        {
            c_bits[idx++] = 0;
        }else if(phy_struct->vd_tb_state[i] > phy_struct->vd_tb_state[i+1]){
            c_bits[idx++] = 1;
        }else{
            // Check to see if the transition has resulted in the same state
            // In this case, if state is 0 then output is 0
            if(phy_struct->vd_tb_state[i] == 0)
            {
                c_bits[idx++] = 0;
            }else{
                c_bits[idx++] = 1;
            }
        }
    }
    *N_c_bits = idx;
}

int main()
{
    PHY_PCFICH *pcfich=new PHY_PCFICH;
    PHY_PHICH *PHICH = new PHY_PHICH;
    FRAME *subframe = new FRAME;
    LIBLTE_PHY *phy_struct = new LIBLTE_PHY;
    RRC_mib_struct *mib= new RRC_mib_struct;
    PDCCH_STRUCT *pdcch=new PDCCH_STRUCT;
    int i;

    //PBCH related
    int in_bits[24]={0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
    int bits[40]={1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1};
    int g[3]={133,177,165};
    int sfn=0;
    int N_id_cell=0;
    int N_ant=4;
    int N_d_bits,N_bits;
    int phich_res;
    bool x;
    x=0;
   // x=mib->phich_config.dur;

    phich_res=mib->phich_config.res;
    init_phy(phy_struct,FS_30_72,N_ant,N_id_cell,PHICH_Ng_value[0],x);
    pbch_channel_map(subframe,phy_struct,&in_bits[0],sfn,N_ant,N_id_cell);
    cout<<"pbch channel mapped"<<endl;

    pcfich_channel_map(phy_struct,N_id_cell,N_ant,pcfich,subframe);
    cout<<"pcfich channel mapped"<<endl;
    phich_channel_map(phy_struct,N_id_cell,N_ant,PHICH,pcfich,subframe);
    cout<<"phich channel mapped"<<endl;
    crs_map(phy_struct,subframe,N_id_cell,N_ant,x);
    cout<<"crs_mapped"<<endl;


   for(i=0;i<phy_struct->N_sc_rb*phy_struct->N_rb_dl;i++)
    {
        cout<<"bch symbol "<<i<<":"<<subframe->tx_symb_re[0][7][i]<<endl;
    }

    for(i=0;i<phy_struct->N_sc_rb*phy_struct->N_rb_dl;i++)
       {
           cout<<"pcfich "<<i<<":"<<subframe->tx_symb_re[0][0][i]<<endl;
       }
    for(i=0;i<phy_struct->N_sc_rb*phy_struct->N_rb_dl;i++)
    {
        cout<<"CRS "<<i<<":"<<subframe->tx_symb_re[0][0][i]<<endl;
    }

    //Dummy : To get a feel of PDCCH allocations
    int N_reg,N_cce;

    N_reg=(pcfich->cfi)*(3*phy_struct->N_rb_dl)-(phy_struct->N_rb_dl);
    if(N_ant==4)
        N_reg-=phy_struct->N_rb_dl;
    N_cce=N_reg/9;

    cout<<"No. of CCE's :"<<N_cce<<endl;
    conv_encode(phy_struct,bits,40,7,3,g,1,&phy_struct->bch_tx_dbits[0],&N_d_bits);
   //viterbi_decode(phy_struct,g,7,3,N_d_bits,&phy_struct->bch_tx_dbits[0],&N_d_bits);
    pdcch_channel_encode(pcfich,PHICH,phy_struct,pdcch,subframe,N_ant,N_id_cell);

}
