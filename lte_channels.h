//
//  phy_channels.h
//  LTE_PHY
//Contains the defines for every channel
//  Created by Vaishnavi  on 08/02/20.
//  Copyright © 2020 none. All rights reserved.
//



#ifndef phy_channels_h
#define phy_channels_h

/*-----------------------------------------------------------------------
 The initialisation function for initially assigning values.

-------------------------------------------------------------------------*/
void init_phy(LIBLTE_PHY *phy_struct,FS_enum fs,int N_ant,int N_id_cell,float phich_res,bool x);
/*-----------------------------------------------------------------------
To calculate CRC for the input data
 -------------------------------------------------------------------------*/
void calc_crc(int  *a_bits,
int  N_a_bits,
int  crc,
int  *p_bits,
int  N_p_bits);
/*-----------------------------------------------------------------------
 Converting the input to bits

 -------------------------------------------------------------------------*/
void convert_to_bits(int val,int **bits,int N_bits);

/*-----------------------------------------------------------------------
 Convolution encoder used by the control channels in the downlink
  -------------------------------------------------------------------------*/
void conv_encode(LIBLTE_PHY *phy_struct,
int             *c_bits,
int             N_c_bits,
int             constraint_len,
int             rate,
int             *g,
bool               tail_bit,
int             *d_bits,
int            *N_d_bits);

/*-----------------------------------------------------------------------
  Used for converting the CFI value into bits for further processing and
 mapping to  pcfich channel.
 Inputs
 N_bits-No. of output bits
 cfi-Value of CFI
 bits-The output bits
--------------------------------------------------------------------------*/
void cfi_encode(LIBLTE_PHY *phy_struct,int cfi,int *bits,int *N_bits);



/*----------------------------------------------------------------------------------
 Used for precoding depending on whether transmit diversity or spatial multiplexing.
 TO FIX:Currently supporting only transmit diversity, to add spatial multiplexing
         Also the precoding for PCHICH for 4 antennas is not correct.(To edit)
 -----------------------------------------------------------------------------------*/
void pre_coder_dl(float        *x_re,
float                          *x_im,
int                          M_layer_symb,
int                           N_ant,
LIBLTE_PHY_PRE_CODER_TYPE_ENUM  type,
float                          *y_re,
float                          *y_im,
int                         y_len,
int                       *M_ap_symb);


/*----------------------------------------------------------------------------------
 Layer mapping is done depending on the no. of antennas.
 Takes the codewords as input and maps to the layers.
 -----------------------------------------------------------------------------------*/
void layer_mapper(float        *d_re,
float                          *d_im,
int                          M_symb,
int                           N_ant,
int                          N_codewords,
LIBLTE_PHY_PRE_CODER_TYPE_ENUM  type,
float                          *x_re,
float                          *x_im,
int                         *M_layer_symb);


/*----------------------------------------------------------------------------------
Modulation mapper is used for modulating the input bits
 depending on the modulation type.
 Takes bits as inputs and converts them into real and imaginary value
-----------------------------------------------------------------------------------*/
void modulation_mapper(int                           *bits,
int                           N_bits,
PHY_MODULATION_TYPE_ENUM  type,
float                           *d_re,
float                           *d_im,
int                          *M_symb);


/*----------------------------------------------------------------------------------
 Used for rate matching for control channels

 -----------------------------------------------------------------------------------*/
void rate_match_conv(LIBLTE_PHY *phy_struct,
int             *d_bits,
int             N_d_bits,
int             N_e_bits,
int             *e_bits);


/*----------------------------------------------------------------------------------
 To generate a scrambling sequence based on the value of c_init.
 The value of c_init depends on the type of control channel
 -----------------------------------------------------------------------------------*/
void scramb_generate(int  c_init,
int len,
int *c);


/*----------------------------------------------------------------------------------
 Generates the CRS sequence for the antenna ports
 ----------------------------------------------------------------------------------*/
void generate_crs(int ns,int l,bool x,int N_id_cell,float *c_re,float *c_im);


/*----------------------------------------------------------------------------------
This maps the generated CRS to the resource grid
----------------------------------------------------------------------------------*/
 void crs_map(LIBLTE_PHY *phy_struct,FRAME *subframe,int N_id_cell,int N_ant,bool x);



/*----------------------------------------------------------------------------------
 Used to convert the input bits at the PHY layer which undergoes PHY layer processing
 and gets mapped to the subframe.This mapping is done every subframe depending on the
 subframe number.
 -----------------------------------------------------------------------------------*/
void pbch_channel_map(FRAME *subframe,LIBLTE_PHY *phy_struct,int *in_bits,int sfn,int N_ant,int N_id_cell);


/*----------------------------------------------------------------------------------
 Used for pcfich channel mapping by taking the input bits,processing and mapping
 to the first symbol of every subframe.
 -----------------------------------------------------------------------------------*/
void pcfich_channel_map(LIBLTE_PHY *phy_struct,int N_id_cell,int N_ant,PHY_PCFICH *pcfich,FRAME *subframe);


/*----------------------------------------------------------------------------------
 Used for phich taking the bits,processing and mapping to the first symbol of every subframe.
 -----------------------------------------------------------------------------------*/
void phich_channel_map(LIBLTE_PHY *phy_struct,int N_id_cell,int N_ant,PHY_PHICH *PHICH,PHY_PCFICH *pcfich,FRAME *subframe);

/*----------------------------------------------------------------------------------
Used for encoding the PDCCH channel for the total no. of allocations,send them through a modulation mapper,
 layer mapper,precoder and mapping them to the resource elements.
 Currently works on aggregation level 4 for both common and user specific search space.
 (Might need to the add the other aggregation levels if needed) Try and find out who decides
 the aggregation level and how does the updation dynamically happen in the PHY layer
-----------------------------------------------------------------------------------*/
void pdcch_channel_encode(PHY_PCFICH *pcfich,PHY_PHICH *phich,LIBLTE_PHY *phy_struct,PDCCH_STRUCT *pdcch,FRAME *subframe,int N_ant,int N_id_cell);

/*----------------------------------------------------------------------------------
 Used for finding the permuted matrix before REG mapping
-----------------------------------------------------------------------------------*/
void pdcch_pre_calc(LIBLTE_PHY *phy_struct,int N_ant);

/*----------------------------------------------------------------------------------
 Encoding the input DCI payload by appending the CRC,XORing with the RNTI and passing
 it to a convolution encoder and rate matcher.
 Note:Only an aggregation level 4 is implemented for both the user and common search
 space for every user.
 FIX ME:The other aggregation levels need to be implemented
 -----------------------------------------------------------------------------------*/
void dci_encode(LIBLTE_PHY *phy_struct,int *in_bits,int rnti,int ue_ant,int N_in_bits,int *out_bits,int N_out_bits);

/*----------------------------------------------------------------------------------
 Used to calculate the no. of CCE's available for PDCCH

 ----------------------------------------------------------------------------------*/
void calculate_cce(PHY_PCFICH *pcfich,PHY_PHICH *phich,LIBLTE_PHY *phy_struct,int N_ant,int pdcch_symbols,int *N_cce,int *N_reg_pdcch);
/*----------------------------------------------------------------------------------
 To generate the DCI payload for DCI0.
 ----------------------------------------------------------------------------------*/
void dci_0_pack(PHY_ALLOCATION_STRUCT *alloc,int N_rb_ul,int N_ant,int *out_bits,int *N_out_bits);

/*----------------------------------------------------------------------------------
To generate DCI payload for DCI 1a
 ----------------------------------------------------------------------------------*/
void dci_1a_pack(PHY_ALLOCATION_STRUCT *alloc,int N_rb_dl,int N_ant,int *out_bits,int *N_out_bits);
#endif /* phy_channels_h */
