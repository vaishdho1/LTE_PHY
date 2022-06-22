//
//  liblte_phy.h
//  LTE_PHY
// Contains all the structures used by the various
//downlink channles to store data.
//  Created by Vaishnavi  on 29/12/19.
//  Copyright © 2019 none. All rights reserved.
//
#include<stdint.h>
#include<fftw3.h>
#include "phy_header.h"
#ifndef liblte_phy_h
#define liblte_phy_h

//Defines
typedef enum{
   DLSCH = 0,
   PCH,
   ULSCH,
   ULCCH,
}PHY_CHAN_TYPE_ENUM;

typedef struct{
PHY_CHAN_TYPE_ENUM       chan_type;
    int                          tbs;
    int                          rv_idx;
    int                          N_prb;
    int                          prb[Nslots_per_subframe][N_rb_max];
    int                          N_codewords;
    int                          N_layers;
    int                          tx_mode;
    int                          harq_retx_count;
    int                          rnti;
    int                           mcs;
    int                           tpc;
    bool                            ndi;
    bool                            dl_alloc;
}PHY_ALLOCATION_STRUCT;

typedef struct{
    PHY_ALLOCATION_STRUCT alloc[N_max_allocations];
    int N_symbs;
    int N_alloc;
}PDCCH_STRUCT;




//Modulation related types
typedef enum{
    MODULATION_TYPE_BPSK = 0,
    MODULATION_TYPE_QPSK,
    MODULATION_TYPE_16QAM,
    MODULATION_TYPE_64QAM,
}PHY_MODULATION_TYPE_ENUM;

typedef struct{

    float rx_symb_re[16][1200];
    float rx_symb_im[16][1200];
    float tx_symb_re[4][16][1200];
    float tx_symb_im[4][16][1200];
    float rx_ce_re[4][16][1200];
    float rx_ce_im[4][16][1200];
    unsigned int num=0;
} FRAME;
//precoder types
typedef enum LIBLTE_PHY_PRE_CODER_TYPE_ENUM{
    spatial=0,
    transmit=1
}LIBLTE_PHY_PRE_CODER_TYPE_ENUM;

typedef struct{
    //Rate match convolution
    int rmc_tmp[MAX_CODE_BLOCK_SIZE];
    int rmc_sb_mat[COLUMNS_RATE_MATCH][MAX_CODE_BLOCK_SIZE/COLUMNS_RATE_MATCH];
    int rmc_sb_perm_mat[COLUMNS_RATE_MATCH][MAX_CODE_BLOCK_SIZE/COLUMNS_RATE_MATCH];
    int rmc_w[BASE_CODING_RATE*MAX_CODE_BLOCK_SIZE];
    int N_rb_dl;
    int N_rb_ul;
    int N_sc_rb=12;
    //bch bits
    int bch_bits[40];
    int bch_N_bits;
    int bch_c[1920];
    int bch_tx_dbits[1920];
    int bch_encode_bits[1920];
    int bch_scramb_bits[1920];
    float bch_x_re[480];
    float bch_x_im[480];
    float bch_d_re[240];
    float bch_d_im[240];
    float bch_y_re[3][240];
    float bch_y_im[3][240];

    //sss related
    float sss_mod_re_0[168][72];
    float sss_mod_im_0[168][72];
    float sss_mod_im_5[168][72];
    float sss_mod_re_5[168][72];
    float sss_re_0[72];
    float sss_re_5[72];
    float sss_im_0[72];
    float sss_im_5[72];

    //pss related
    float pss_mod_re[3][72];
    float pss_mod_im[3][72];
    float pss_mod_re_n1[3][72];
    float pss_mod_im_n1[3][72];
    float pss_mod_re_p1[3][72];
    float pss_mod_im_p1[3][72];
    float rx_symb_re[72];
    float rx_symb_im[72];

    //CRS related
    float crs_re_storage[20][3][N_rb_max*N_sc_rb_max];
    float crs_im_storage[20][3][N_rb_max*N_sc_rb_max];
    float crs_re[2*N_symbols_max][N_rb_max*N_sc_rb_max];
    float crs_im[2*N_symbols_max][N_rb_max*N_sc_rb_max];
    int N_id_cell_crs;


    //PHICH

    int N_sf_phich;//Depends on the cyclic prefix length-normal:4,extended:2
    int N_group_phich;//Depends on the value of Ng and N_rb_dl
    fftw_complex *s2s_in;
    fftw_complex *s2s_out;
    fftw_plan symbs_to_samps_dl_plan;
    fftw_plan samps_to_symbs_dl_plan;
    //PDCCH
    float pdcch_re[3][N_CCE_MAX][576];
    float pdcch_im[3][N_CCE_MAX][576];
    float pdcch_reg_re[3][N_REG_MAX][576];
    float pdcch_reg_im[3][N_REG_MAX][4];
    float pdcch_x_re[576];
    float pdcch_x_im[576];
    float pdcch_y_re[3][576];
    float pdcch_y_im[3][576];
    float pdcch_d_re[576];
    float pdcch_d_im[576];
    float pdcch_perm_re[3][N_REG_MAX][4];
    float pdcch_perm_im[3][N_REG_MAX][4];
    float pdcch_shift_re[3][N_REG_MAX][4];
    float pdcch_shift_im[3][N_REG_MAX][4];
    int pdcch_encode_bits[576];
    int pdcch_scramb_bits[576];
    int pdcch_c[576*2];
    int pdcch_reg_vec[N_REG_MAX];
    int pdcch_reg_perm_vec[N_REG_MAX];
    int pdcch_permute_map[N_REG_MAX][N_REG_MAX];
    int pdcch_cce_re[N_CCE_MAX][N_REG_CCE*4];
    int pdcch_cce_im[N_CCE_MAX][N_REG_CCE*4];
    bool pdcch_used[N_CCE_MAX];
    int pdcch_y_re_est[PDCCH_BITS_MAX];
    int pdcch_y_im_est[PDCCH_BITS_MAX];
    int pdcch_c_re_est[4][PDCCH_BITS_MAX/2];
    int pdcch_c_im_est[4][PDCCH_BITS_MAX/2];
    int pdcch_dci[PDCCH_BITS_MAX];
    //DCI

    int dci_rx_d_bits[576];
    int dci_tx_d_bits[576];
    int dci_c[192];


    //General configs
    int N_ant;
    int N_id_cell;
    int N_samps_per_slot;
    int N_samps_per_frame;
    int N_samps_cp_1_else;
    int N_samps_cp_1_0;
    int N_samps_per_symb;
    int N_samps_per_subfr;
    int fs;
    //Channel estimation

    float dl_ce_crs_re[16][N_rb_max*N_sc_rb_max];
    float dl_ce_crs_im[16][N_rb_max*N_sc_rb_max];
    float dl_ce_mag[16][N_rb_max*N_sc_rb_max];
    float dl_ce_ang[16][N_rb_max*N_sc_rb_max];


    //Viterbi decoding
    float vd_path_metric[LIBLTE_PHY_MAX_VITERBI_STATES][MAX_CODE_BLOCK_SIZE];
    float vd_br_metric[LIBLTE_PHY_MAX_VITERBI_STATES][2];
    float vd_p_metric[LIBLTE_PHY_MAX_VITERBI_STATES][2];
    float vd_br_weight[LIBLTE_PHY_MAX_VITERBI_STATES][2];
    float vd_w_metric[LIBLTE_PHY_MAX_VITERBI_STATES][MAX_CODE_BLOCK_SIZE];
    float vd_tb_state[MAX_CODE_BLOCK_SIZE];
    float vd_tb_weight[MAX_CODE_BLOCK_SIZE];
    int vd_st_output[LIBLTE_PHY_MAX_VITERBI_STATES][2][3];

}LIBLTE_PHY;

typedef struct{

    int N_reg;

    bool present[2][8]={{1,1,0,0,1,0,1,1},{1,1,1,1,0,0,0,1}};//To be determined
    int b[2][8]={{1,1,1,1,1,0,1,0},{1,1,1,0,0,0,0,1}}; //Obtained from the higher layers
    int k[4];
    float z_re[12];
    float z_im[12];

} PHY_PHICH;

typedef struct{
    float n[4];
    int k[4];
    int N_reg=4; //Constant for pcfich;
    int cfi = 3; //Need to find where exactly it is found(Looks like programmed by the user)

}PHY_PCFICH;

typedef enum{
    FS_1_92=0,
    FS_3_84,
    FS_7_68,
    FS_15_36,
    FS_30_72,

}FS_enum;

#endif /* liblte_phy_h */
