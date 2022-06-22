//
//  RRC_header.h
//  LTE_PHY
// Contains RRC related structures
//  Created by Vaishnavi  on 07/02/20.
//  Copyright © 2020 none. All rights reserved.
//

#ifndef RRC_header_h
#define RRC_header_h


//MIB information
typedef enum{
    PHICH_normal=0,
    PHICH_extended,
}PHICH_duration_enum;

typedef enum{
    PHICH_Ng_1_6=0,
    PHICH_Ng_1_2,
    PHICH_Ng_1,
    PHICH_Ng_2,
    PHICH_N_items,

}PHICH_res_enum;

typedef struct{
    PHICH_duration_enum dur;
    PHICH_res_enum res;
}RRC_PHICH_config;

    typedef enum{
        DL_BW_1_4=0,
        DL_BW_6,
        DL_BW_15,
        DL_BW_25,
        DL_BW_50,
        DL_BW_75,
        DL_BW_100,
    }RRC_DL_BW_enum;


typedef struct{
    RRC_PHICH_config phich_config;
    RRC_DL_BW_enum dl_bw;
    int sfn_div_4;//This is the divided by 4 version of sfn as only 8bits are used to encode this in the system

}RRC_mib_struct;






static double PHICH_Ng_value[PHICH_N_items]={.1666667,.5,1,2};




#endif /* RRC_header_h */
