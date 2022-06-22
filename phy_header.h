//
//  phy_header.h
//  LTE_PHY
//  Macros defined for various channel related configurations
//  Created by Vaishnavi  on 07/02/20.
//  Copyright © 2020 none. All rights reserved.
//

#ifndef phy_header_h
#define phy_header_h
//-------------DEFINES---------------------//

#define DCI_0_1A_FLAG_0          0
#define DCI_0_1A_FLAG_1A         1
#define DCI_VRB_TYPE_LOCALIZED   0
#define DCI_VRB_TYPE_DISTRIBUTED 1


//CRC
#define CRC16  0x00011021
#define CRC24A 0x01864CFB
#define CRC24B 0x01800063
#define CRC8   0x0000019B

#define MAX_CODE_BLOCK_SIZE 6176
#define COLUMNS_RATE_MATCH 32
#define BASE_CODING_RATE 3

#define TX_NULL_SYMB 100

//PDCCH DEFINES
#define N_max_allocations 6
#define N_REG_MAX 767
#define N_CCE_MAX N_REG_MAX/N_REG_CCE
#define N_REG_CCE 9
#define PDCCH_BITS_MAX 576 //Corresponding to PDCCH format 3

//RNTI values

#define MAC_SI_RNTI_START 0x01
#define MAC_SI_RNTI_END 0x3C
#define MAC_CRNTI_START 0x3D
#define MAC_CRNTI_END 0xFFF3
#define MAC_RESV_RNTI_START 0xFFF4
#define MAC_RESV_RNTI_END 0xFFFC
#define MAC_MRNTI 0xFFFD
#define MAC_PRNTI 0xFFFE
#define MAC_SIRNTI 0xFFFF
//Viterbi algorithm


#define LIBLTE_PHY_MAX_VITERBI_STATES 128


//CFI (PCFICH related)

int CFI_0[32]={0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,0,1};
int CFI_1[32]={1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0};
int CFI_2[32]={1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1};
int CFI_3[32]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

//CRC calculation

int ICC_PERM[32]={1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31,
0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30};

float PHICH_normal_re[8][4]={{1,1,1,1},{1,1,1,-1},{1,1,-1,-1},{1,-1,-1,1},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
float PHICH_normal_im[8][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{1,1,1,1},{1,-1,1,-1},{1,1,-1,-1},{1,-1,-1,1}};

float PHICH_extended_re[4][2]={{1,1},{1,-1},{0,0},{0,0}};
float PHICH_extended_im[4][2]={{0,0},{0,0},{1,1},{1,-1}};


#define Nslots_per_subframe 2
#define N_subframe_per_frame 10
#define N_sc_rb_max 12
#define N_rb_max 100
#define N_symbols_max 7

//Bandwidth and no. of samples
//1.4 MHz


#define N_samps_per_symbol_1_92 128
#define N_samples_slot_1_92 960
#define N_samps_CP_1_92 10
#define N_samps_CP_else_1_92 9
#define N_samps_per_subframe_1_92 (N_samples_slot_1_92*Nslots_per_subframe)
#define N_samps_per_frame_1_92 (N_samps_per_subframe_1_92*N_subframe_per_frame)
//3 MHz
#define N_samps_per_symbol_3_84 256
#define N_samples_slot_3_84 1920
#define N_samps_CP_3_84 20
#define N_samps_CP_else_3_84 18
#define N_samps_per_subframe_3_84 (N_samples_slot_3_84*Nslots_per_subframe)
#define N_samps_per_frame_3_84 (N_samps_per_subframe_3_84*N_subframe_per_frame)
//5 MHz
#define N_samps_per_symbol_7_68 512
#define N_samples_slot_7_68 3840
#define N_samps_CP_7_68 40
#define N_samps_CP_else_7_68 36
#define N_samps_per_subframe_7_68 (N_samples_slot_7_68*Nslots_per_subframe)
#define N_samps_per_frame_7_68 (N_samps_per_subframe_7_68*N_subframe_per_frame)
//10 MHz
#define N_samps_per_symbol_15_36 1024
#define N_samples_slot_15_36 7680
#define N_samps_CP_15_36 80
#define N_samps_CP_else_15_36 72
#define N_samps_per_subframe_15_36 (N_samples_slot_15_36*Nslots_per_subframe)
#define N_samps_per_frame_15_36 (N_samps_per_subframe_15_36*N_subframe_per_frame)
//15 & 20MHz
#define N_samps_per_symbol_30_72 2048
#define N_samples_slot_30_72 15360
#define N_samps_CP_30_72 160
#define N_samps_CP_else_30_72 144
#define N_samps_per_subframe_30_72 (N_samples_slot_30_72*Nslots_per_subframe)
#define N_samps_per_frame_30_72 (N_samps_per_subframe_30_72*N_subframe_per_frame)


#define N_REG_CCE 9



#endif /* phy_header_h */
