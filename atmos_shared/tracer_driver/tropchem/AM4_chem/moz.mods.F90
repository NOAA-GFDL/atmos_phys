      module mo_grid_mod
!---------------------------------------------------------------------
! ... Basic grid point resolution parameters
!---------------------------------------------------------------------
      implicit none
      save
      integer, parameter :: &
                pcnst = 122 +1, & ! number of advected constituents including cloud water
                pcnstm1 = 122, & ! number of advected constituents excluding cloud water
                plev = 1, & ! number of vertical levels
                plevp = plev+1, & ! plev plus 1
                plevm = plev-1, & ! plev minus 1
                plon = 1, & ! number of longitudes
                plat = 1 ! number of latitudes
      integer, parameter :: &
                pnats = 0 ! number of non-advected trace species
      integer :: nodes ! mpi task count
      integer :: plonl ! longitude tile dimension
      integer :: pplon ! longitude tile count
      integer :: plnplv ! plonl * plev
      end module mo_grid_mod
      module chem_mods_mod
!--------------------------------------------------------------
! ... basic chemistry array parameters
!--------------------------------------------------------------
      use mo_grid_mod, only : pcnstm1
      use mpp_mod, only : mpp_error, FATAL
      implicit none
      save
      integer, parameter :: hetcnt = 0, & ! number of heterogeneous processes
                            phtcnt = 48, & ! number of photo processes
                            rxntot = 305, & ! number of total reactions
                            gascnt = 257, & ! number of gas phase reactions
                            nfs = 3, & ! number of "fixed" species
                            relcnt = 0, & ! number of relationship species
                            grpcnt = 0, & ! number of group members
                            imp_nzcnt = 1360, & ! number of non-zero implicit matrix entries
                            rod_nzcnt = 0, & ! number of non-zero rodas matrix entries
                            extcnt = 0, & ! number of species with external forcing
                            clscnt1 = 0, & ! number of species in explicit class
                            clscnt2 = 0, & ! number of species in hov class
                            clscnt3 = 0, & ! number of species in ebi class
                            clscnt4 = 122, & ! number of species in implicit class
                            clscnt5 = 0, & ! number of species in rodas class
                            indexm = 1, & ! index of total atm density in invariant array
                            ncol_abs = 2, & ! number of column densities
                            indexh2o = 0, & ! index of water vapor density
                            clsze = 1 ! loop length for implicit chemistry
      integer :: ngrp = 0
      integer :: drydep_cnt = 0
      integer :: srfems_cnt = 0
      integer :: rxt_alias_cnt = 0
      integer, allocatable :: grp_mem_cnt(:)
      integer, allocatable :: rxt_alias_map(:)
      real :: adv_mass(pcnstm1)
      real :: nadv_mass(grpcnt)
      character(len=16), allocatable :: rxt_alias_lst(:)
      character(len=8), allocatable :: drydep_lst(:)
      character(len=8), allocatable :: srfems_lst(:)
      character(len=8), allocatable :: grp_lst(:)
      character(len=8) :: het_lst(max(1,hetcnt))
      character(len=8) :: extfrc_lst(max(1,extcnt))
      character(len=8) :: inv_lst(max(1,nfs))
      type solver_class
  integer :: clscnt
  integer :: lin_rxt_cnt
  integer :: nln_rxt_cnt
  integer :: indprd_cnt
  integer :: iter_max
         integer :: cls_rxt_cnt(4)
         integer, pointer :: permute(:)
         integer, pointer :: diag_map(:)
         integer, pointer :: clsmap(:)
      end type solver_class
      type(solver_class) :: explicit, implicit, rodas
      contains
      subroutine endrun(msg)
      implicit none
      character(len=128), intent(in), optional :: msg
      call mpp_error(FATAL, msg)
      end subroutine endrun
      subroutine chem_mods_init
!--------------------------------------------------------------
! ... intialize the class derived type
!--------------------------------------------------------------
      implicit none
      integer :: astat
      explicit%clscnt = 0
      explicit%indprd_cnt = 0
      implicit%clscnt = 122
      implicit%lin_rxt_cnt = 96
      implicit%nln_rxt_cnt = 208
      implicit%indprd_cnt = 1
      implicit%iter_max = 11
      rodas%clscnt = 0
      rodas%lin_rxt_cnt = 0
      rodas%nln_rxt_cnt = 0
      rodas%indprd_cnt = 0
      if( explicit%clscnt > 0 ) then
  allocate( explicit%clsmap(explicit%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate explicit%clsmap ; error = ',astat
     call endrun
  end if
         explicit%clsmap(:) = 0
      end if
      if( implicit%clscnt > 0 ) then
  allocate( implicit%permute(implicit%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate implicit%permute ; error = ',astat
     call endrun
  end if
         implicit%permute(:) = 0
  allocate( implicit%diag_map(implicit%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate implicit%diag_map ; error = ',astat
     call endrun
  end if
         implicit%diag_map(:) = 0
  allocate( implicit%clsmap(implicit%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate implicit%clsmap ; error = ',astat
     call endrun
  end if
         implicit%clsmap(:) = 0
      end if
      if( rodas%clscnt > 0 ) then
  allocate( rodas%permute(rodas%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate rodas%permute ; error = ',astat
     call endrun
  end if
         rodas%permute(:) = 0
  allocate( rodas%diag_map(rodas%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate rodas%diag_map ; error = ',astat
     call endrun
  end if
         rodas%diag_map(:) = 0
  allocate( rodas%clsmap(rodas%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate rodas%clsmap ; error = ',astat
     call endrun
  end if
         rodas%clsmap(:) = 0
      end if
      end subroutine chem_mods_init
      end module chem_mods_mod
      module M_SPC_ID_MOD
      implicit none
      integer, parameter :: id_O3 = 1
      integer, parameter :: id_O = 2
      integer, parameter :: id_O1D = 3
      integer, parameter :: id_N2O = 4
      integer, parameter :: id_N = 5
      integer, parameter :: id_NO = 6
      integer, parameter :: id_NO2 = 7
      integer, parameter :: id_NO3 = 8
      integer, parameter :: id_HNO3 = 9
      integer, parameter :: id_HO2NO2 = 10
      integer, parameter :: id_N2O5 = 11
      integer, parameter :: id_CH4 = 12
      integer, parameter :: id_CH3O2 = 13
      integer, parameter :: id_HNO3_D1 = 14
      integer, parameter :: id_HNO3_D2 = 15
      integer, parameter :: id_HNO3_D3 = 16
      integer, parameter :: id_HNO3_D4 = 17
      integer, parameter :: id_HNO3_D5 = 18
      integer, parameter :: id_SO4_D1 = 19
      integer, parameter :: id_SO4_D2 = 20
      integer, parameter :: id_SO4_D3 = 21
      integer, parameter :: id_SO4_D4 = 22
      integer, parameter :: id_SO4_D5 = 23
      integer, parameter :: id_CH3OOH = 24
      integer, parameter :: id_CH2O = 25
      integer, parameter :: id_CO = 26
      integer, parameter :: id_OH = 27
      integer, parameter :: id_HO2 = 28
      integer, parameter :: id_H2O2 = 29
      integer, parameter :: id_C3H6 = 30
      integer, parameter :: id_ISOP = 31
      integer, parameter :: id_PO2 = 32
      integer, parameter :: id_CH3CHO = 33
      integer, parameter :: id_POOH = 34
      integer, parameter :: id_CH3CO3 = 35
      integer, parameter :: id_CH3COOOH = 36
      integer, parameter :: id_PAN = 37
      integer, parameter :: id_C2H6 = 38
      integer, parameter :: id_C2H4 = 39
      integer, parameter :: id_C4H10 = 40
      integer, parameter :: id_MPAN = 41
      integer, parameter :: id_ISOPO2 = 42
      integer, parameter :: id_MVK = 43
      integer, parameter :: id_MACR = 44
      integer, parameter :: id_MACRO2 = 45
      integer, parameter :: id_MACROOH = 46
      integer, parameter :: id_C2H5O2 = 47
      integer, parameter :: id_C2H5OOH = 48
      integer, parameter :: id_C10H16 = 49
      integer, parameter :: id_C3H8 = 50
      integer, parameter :: id_C3H7O2 = 51
      integer, parameter :: id_C3H7OOH = 52
      integer, parameter :: id_CH3COCH3 = 53
      integer, parameter :: id_CH3OH = 54
      integer, parameter :: id_C2H5OH = 55
      integer, parameter :: id_GLYALD = 56
      integer, parameter :: id_HYAC = 57
      integer, parameter :: id_EO2 = 58
      integer, parameter :: id_EO = 59
      integer, parameter :: id_ISOPOOH = 60
      integer, parameter :: id_H2 = 61
      integer, parameter :: id_SO2 = 62
      integer, parameter :: id_SO4 = 63
      integer, parameter :: id_DMS = 64
      integer, parameter :: id_NH3 = 65
      integer, parameter :: id_NH4NO3 = 66
      integer, parameter :: id_NH4 = 67
      integer, parameter :: id_HCl = 68
      integer, parameter :: id_HOCl = 69
      integer, parameter :: id_ClONO2 = 70
      integer, parameter :: id_Cl = 71
      integer, parameter :: id_ClO = 72
      integer, parameter :: id_Cl2O2 = 73
      integer, parameter :: id_Cl2 = 74
      integer, parameter :: id_HOBr = 75
      integer, parameter :: id_HBr = 76
      integer, parameter :: id_BrONO2 = 77
      integer, parameter :: id_Br = 78
      integer, parameter :: id_BrO = 79
      integer, parameter :: id_BrCl = 80
      integer, parameter :: id_LCH4 = 81
      integer, parameter :: id_H = 82
      integer, parameter :: id_H2O = 83
      integer, parameter :: id_ROH = 84
      integer, parameter :: id_RCHO = 85
      integer, parameter :: id_RCO3 = 86
      integer, parameter :: id_ISOPNB = 87
      integer, parameter :: id_ISOPNBO2 = 88
      integer, parameter :: id_MACRN = 89
      integer, parameter :: id_MVKN = 90
      integer, parameter :: id_R4N2 = 91
      integer, parameter :: id_MEK = 92
      integer, parameter :: id_MEKO2 = 93
      integer, parameter :: id_MEKOOH = 94
      integer, parameter :: id_R4N1 = 95
      integer, parameter :: id_IEPOX = 96
      integer, parameter :: id_IEPOXOO = 97
      integer, parameter :: id_GLYX = 98
      integer, parameter :: id_MGLY = 99
      integer, parameter :: id_MVKO2 = 100
      integer, parameter :: id_MVKOOH = 101
      integer, parameter :: id_MACRNO2 = 102
      integer, parameter :: id_MAO3 = 103
      integer, parameter :: id_MAOP = 104
      integer, parameter :: id_MAOPO2 = 105
      integer, parameter :: id_ATO2 = 106
      integer, parameter :: id_ATOOH = 107
      integer, parameter :: id_INO2 = 108
      integer, parameter :: id_INPN = 109
      integer, parameter :: id_ISNOOA = 110
      integer, parameter :: id_ISN1 = 111
      integer, parameter :: id_O3S = 112
      integer, parameter :: id_O3S_E90 = 113
      integer, parameter :: id_CH3SH = 114
      integer, parameter :: id_CH3S = 115
      integer, parameter :: id_CH3SO = 116
      integer, parameter :: id_CH3SOO = 117
      integer, parameter :: id_CH3SO2 = 118
      integer, parameter :: id_CH3SO3 = 119
      integer, parameter :: id_MSA = 120
      integer, parameter :: id_MTMP = 121
      integer, parameter :: id_HPMTF = 122
      end module M_SPC_ID_MOD
      module M_RXT_ID_MOD
      implicit none
      integer, parameter :: rid_jo2 = 1
      integer, parameter :: rid_jo1d = 2
      integer, parameter :: rid_jo3p = 3
      integer, parameter :: rid_jn2o = 4
      integer, parameter :: rid_jno = 5
      integer, parameter :: rid_jno2 = 6
      integer, parameter :: rid_jn2o5 = 7
      integer, parameter :: rid_jhno3 = 8
      integer, parameter :: rid_jno3 = 9
      integer, parameter :: rid_jho2no2 = 10
      integer, parameter :: rid_jch3ooh = 11
      integer, parameter :: rid_jch2o_a = 12
      integer, parameter :: rid_jch2o_b = 13
      integer, parameter :: rid_jh2o = 14
      integer, parameter :: rid_jh2o2 = 15
      integer, parameter :: rid_jch3cho = 16
      integer, parameter :: rid_jpooh = 17
      integer, parameter :: rid_jch3co3h = 18
      integer, parameter :: rid_jpan = 19
      integer, parameter :: rid_jmpan = 20
      integer, parameter :: rid_jmacr_a = 21
      integer, parameter :: rid_jmacr_b = 22
      integer, parameter :: rid_jmvk = 23
      integer, parameter :: rid_jc2h5ooh = 24
      integer, parameter :: rid_jc3h7ooh = 25
      integer, parameter :: rid_jacet = 26
      integer, parameter :: rid_jmgly = 27
      integer, parameter :: rid_jglyoxal1 = 28
      integer, parameter :: rid_jglyoxal2 = 29
      integer, parameter :: rid_jglyoxal3 = 30
      integer, parameter :: rid_jisopooh = 31
      integer, parameter :: rid_jhyac = 32
      integer, parameter :: rid_jglyald = 33
      integer, parameter :: rid_jisopnb = 34
      integer, parameter :: rid_jmacrn = 35
      integer, parameter :: rid_jmvkn = 36
      integer, parameter :: rid_jr4n2 = 37
      integer, parameter :: rid_jmek = 38
      integer, parameter :: rid_jmekooh = 39
      integer, parameter :: rid_jrcho = 40
      integer, parameter :: rid_jclono2 = 41
      integer, parameter :: rid_jhocl = 42
      integer, parameter :: rid_jcl2o2 = 43
      integer, parameter :: rid_jbrono2 = 44
      integer, parameter :: rid_jhobr = 45
      integer, parameter :: rid_jbrcl = 46
      integer, parameter :: rid_jbro = 47
      integer, parameter :: rid_jcl2 = 48
      integer, parameter :: rid_uo_o2 = 49
      integer, parameter :: rid_uco_oha = 53
      integer, parameter :: rid_uco_ohb = 54
      integer, parameter :: rid_ol_oh = 58
      integer, parameter :: rid_ol_ho2 = 59
      integer, parameter :: rid_uho2_ho2 = 60
      integer, parameter :: rid_o1d_n2 = 65
      integer, parameter :: rid_o1d_o2 = 66
      integer, parameter :: rid_ol_o1d = 67
      integer, parameter :: rid_op_ho2 = 70
      integer, parameter :: rid_uno2_no3 = 75
      integer, parameter :: rid_un2o5 = 76
      integer, parameter :: rid_uoh_no2 = 78
      integer, parameter :: rid_uoh_hno3 = 79
      integer, parameter :: rid_uho2_no2 = 81
      integer, parameter :: rid_uhno4 = 83
      integer, parameter :: rid_op_mo2 = 86
      integer, parameter :: rid_uoh_c2h4 = 93
      integer, parameter :: rid_op_eo2 = 94
      integer, parameter :: rid_ol_c2h4 = 97
      integer, parameter :: rid_uoh_c3h6 = 98
      integer, parameter :: rid_ol_c3h6 = 99
      integer, parameter :: rid_op_po2 = 101
      integer, parameter :: rid_op_ch3co3 = 106
      integer, parameter :: rid_upan_f = 107
      integer, parameter :: rid_upan_b = 111
      integer, parameter :: rid_op_c2h5o2 = 114
      integer, parameter :: rid_op_c3h7o2 = 121
      integer, parameter :: rid_uoh_acet = 125
      integer, parameter :: rid_ol_isop = 130
      integer, parameter :: rid_op_isopo2 = 131
      integer, parameter :: rid_op_isopnbo2 = 136
      integer, parameter :: rid_ol_isopnb = 138
      integer, parameter :: rid_op_iepoxo2 = 143
      integer, parameter :: rid_ol_mvk = 145
      integer, parameter :: rid_op_mvko2 = 146
      integer, parameter :: rid_ol_macr = 153
      integer, parameter :: rid_op_macro2 = 155
      integer, parameter :: rid_op_macrno2 = 161
      integer, parameter :: rid_op_mao3 = 163
      integer, parameter :: rid_umpan_f = 167
      integer, parameter :: rid_umpan_b = 168
      integer, parameter :: rid_op_maopo2 = 173
      integer, parameter :: rid_op_ato2 = 175
      integer, parameter :: rid_op_ino2 = 187
      integer, parameter :: rid_ol_c10h16 = 200
      integer, parameter :: rid_uoh_mek = 202
      integer, parameter :: rid_op_meko2 = 203
      integer, parameter :: rid_uoh_rcho = 206
      integer, parameter :: rid_op_rco3 = 208
      integer, parameter :: rid_n2o5h = 225
      integer, parameter :: rid_no3h = 226
      integer, parameter :: rid_ho2h = 227
      integer, parameter :: rid_no2h = 228
      integer, parameter :: rid_so2h = 229
      integer, parameter :: rid_uoh_dms = 232
      integer, parameter :: rid_mtmp_i = 239
      integer, parameter :: rid_hpmtf_a = 240
      integer, parameter :: rid_hno3_d1 = 242
      integer, parameter :: rid_hno3_d2 = 243
      integer, parameter :: rid_hno3_d3 = 244
      integer, parameter :: rid_hno3_d4 = 245
      integer, parameter :: rid_hno3_d5 = 246
      integer, parameter :: rid_no3_d1 = 247
      integer, parameter :: rid_no3_d2 = 248
      integer, parameter :: rid_no3_d3 = 249
      integer, parameter :: rid_no3_d4 = 250
      integer, parameter :: rid_no3_d5 = 251
      integer, parameter :: rid_n2o5_d1 = 252
      integer, parameter :: rid_n2o5_d2 = 253
      integer, parameter :: rid_n2o5_d3 = 254
      integer, parameter :: rid_n2o5_d4 = 255
      integer, parameter :: rid_n2o5_d5 = 256
      integer, parameter :: rid_so2_d1 = 257
      integer, parameter :: rid_so2_d2 = 258
      integer, parameter :: rid_so2_d3 = 259
      integer, parameter :: rid_so2_d4 = 260
      integer, parameter :: rid_so2_d5 = 261
      integer, parameter :: rid_nh3h = 262
      integer, parameter :: rid_strat13 = 264
      integer, parameter :: rid_strat14 = 265
      integer, parameter :: rid_strat20 = 266
      integer, parameter :: rid_strat21 = 267
      integer, parameter :: rid_strat22 = 268
      integer, parameter :: rid_strat23 = 269
      integer, parameter :: rid_strat24 = 270
      integer, parameter :: rid_strat25 = 271
      integer, parameter :: rid_strat26 = 272
      integer, parameter :: rid_strat27 = 273
      integer, parameter :: rid_strat28 = 274
      integer, parameter :: rid_strat29 = 275
      integer, parameter :: rid_strat33 = 276
      integer, parameter :: rid_strat35 = 277
      integer, parameter :: rid_strat37 = 278
      integer, parameter :: rid_strat38 = 279
      integer, parameter :: rid_strat39 = 280
      integer, parameter :: rid_strat40 = 281
      integer, parameter :: rid_strat41 = 282
      integer, parameter :: rid_strat42 = 283
      integer, parameter :: rid_strat43 = 284
      integer, parameter :: rid_strat44 = 285
      integer, parameter :: rid_strat45 = 286
      integer, parameter :: rid_strat46 = 287
      integer, parameter :: rid_strat47 = 288
      integer, parameter :: rid_strat48 = 289
      integer, parameter :: rid_strat69 = 290
      integer, parameter :: rid_strat58 = 291
      integer, parameter :: rid_strat59 = 292
      integer, parameter :: rid_strat64 = 293
      integer, parameter :: rid_strat71 = 294
      integer, parameter :: rid_strat72 = 295
      integer, parameter :: rid_strat73 = 296
      integer, parameter :: rid_strat74 = 297
      integer, parameter :: rid_strat75 = 298
      integer, parameter :: rid_strat76 = 299
      integer, parameter :: rid_strat77 = 300
      integer, parameter :: rid_strat78 = 301
      integer, parameter :: rid_strat79 = 302
      integer, parameter :: rid_strat80 = 303
      integer, parameter :: rid_r0050 = 50
      integer, parameter :: rid_r0051 = 51
      integer, parameter :: rid_r0052 = 52
      integer, parameter :: rid_r0055 = 55
      integer, parameter :: rid_r0056 = 56
      integer, parameter :: rid_r0057 = 57
      integer, parameter :: rid_r0061 = 61
      integer, parameter :: rid_r0062 = 62
      integer, parameter :: rid_r0063 = 63
      integer, parameter :: rid_r0064 = 64
      integer, parameter :: rid_r0068 = 68
      integer, parameter :: rid_r0069 = 69
      integer, parameter :: rid_r0071 = 71
      integer, parameter :: rid_r0072 = 72
      integer, parameter :: rid_r0073 = 73
      integer, parameter :: rid_r0074 = 74
      integer, parameter :: rid_r0077 = 77
      integer, parameter :: rid_r0080 = 80
      integer, parameter :: rid_r0082 = 82
      integer, parameter :: rid_r0084 = 84
      integer, parameter :: rid_r0085 = 85
      integer, parameter :: rid_r0087 = 87
      integer, parameter :: rid_r0088 = 88
      integer, parameter :: rid_r0089 = 89
      integer, parameter :: rid_r0090 = 90
      integer, parameter :: rid_r0091 = 91
      integer, parameter :: rid_r0092 = 92
      integer, parameter :: rid_r0095 = 95
      integer, parameter :: rid_r0096 = 96
      integer, parameter :: rid_r0100 = 100
      integer, parameter :: rid_r0102 = 102
      integer, parameter :: rid_r0103 = 103
      integer, parameter :: rid_r0104 = 104
      integer, parameter :: rid_r0105 = 105
      integer, parameter :: rid_r0108 = 108
      integer, parameter :: rid_r0109 = 109
      integer, parameter :: rid_r0110 = 110
      integer, parameter :: rid_r0112 = 112
      integer, parameter :: rid_r0113 = 113
      integer, parameter :: rid_r0115 = 115
      integer, parameter :: rid_r0116 = 116
      integer, parameter :: rid_r0117 = 117
      integer, parameter :: rid_r0118 = 118
      integer, parameter :: rid_r0119 = 119
      integer, parameter :: rid_r0120 = 120
      integer, parameter :: rid_r0122 = 122
      integer, parameter :: rid_r0123 = 123
      integer, parameter :: rid_r0124 = 124
      integer, parameter :: rid_r0126 = 126
      integer, parameter :: rid_r0127 = 127
      integer, parameter :: rid_r0128 = 128
      integer, parameter :: rid_r0129 = 129
      integer, parameter :: rid_r0132 = 132
      integer, parameter :: rid_r0133 = 133
      integer, parameter :: rid_r0134 = 134
      integer, parameter :: rid_r0135 = 135
      integer, parameter :: rid_r0137 = 137
      integer, parameter :: rid_r0139 = 139
      integer, parameter :: rid_r0140 = 140
      integer, parameter :: rid_r0141 = 141
      integer, parameter :: rid_r0142 = 142
      integer, parameter :: rid_r0144 = 144
      integer, parameter :: rid_r0147 = 147
      integer, parameter :: rid_r0148 = 148
      integer, parameter :: rid_r0149 = 149
      integer, parameter :: rid_r0150 = 150
      integer, parameter :: rid_r0151 = 151
      integer, parameter :: rid_r0152 = 152
      integer, parameter :: rid_r0154 = 154
      integer, parameter :: rid_r0156 = 156
      integer, parameter :: rid_r0157 = 157
      integer, parameter :: rid_r0158 = 158
      integer, parameter :: rid_r0159 = 159
      integer, parameter :: rid_r0160 = 160
      integer, parameter :: rid_r0162 = 162
      integer, parameter :: rid_r0164 = 164
      integer, parameter :: rid_r0165 = 165
      integer, parameter :: rid_r0166 = 166
      integer, parameter :: rid_r0169 = 169
      integer, parameter :: rid_r0170 = 170
      integer, parameter :: rid_r0171 = 171
      integer, parameter :: rid_r0172 = 172
      integer, parameter :: rid_r0174 = 174
      integer, parameter :: rid_r0176 = 176
      integer, parameter :: rid_r0177 = 177
      integer, parameter :: rid_r0178 = 178
      integer, parameter :: rid_r0179 = 179
      integer, parameter :: rid_r0180 = 180
      integer, parameter :: rid_r0181 = 181
      integer, parameter :: rid_r0182 = 182
      integer, parameter :: rid_r0183 = 183
      integer, parameter :: rid_r0184 = 184
      integer, parameter :: rid_r0185 = 185
      integer, parameter :: rid_r0186 = 186
      integer, parameter :: rid_r0188 = 188
      integer, parameter :: rid_r0189 = 189
      integer, parameter :: rid_r0190 = 190
      integer, parameter :: rid_r0191 = 191
      integer, parameter :: rid_r0192 = 192
      integer, parameter :: rid_r0193 = 193
      integer, parameter :: rid_r0194 = 194
      integer, parameter :: rid_r0195 = 195
      integer, parameter :: rid_r0196 = 196
      integer, parameter :: rid_r0197 = 197
      integer, parameter :: rid_r0198 = 198
      integer, parameter :: rid_r0199 = 199
      integer, parameter :: rid_r0201 = 201
      integer, parameter :: rid_r0204 = 204
      integer, parameter :: rid_r0205 = 205
      integer, parameter :: rid_r0207 = 207
      integer, parameter :: rid_r0209 = 209
      integer, parameter :: rid_r0210 = 210
      integer, parameter :: rid_r0211 = 211
      integer, parameter :: rid_r0212 = 212
      integer, parameter :: rid_r0213 = 213
      integer, parameter :: rid_r0214 = 214
      integer, parameter :: rid_r0215 = 215
      integer, parameter :: rid_r0216 = 216
      integer, parameter :: rid_r0217 = 217
      integer, parameter :: rid_r0218 = 218
      integer, parameter :: rid_r0219 = 219
      integer, parameter :: rid_r0220 = 220
      integer, parameter :: rid_r0221 = 221
      integer, parameter :: rid_r0222 = 222
      integer, parameter :: rid_r0223 = 223
      integer, parameter :: rid_r0224 = 224
      integer, parameter :: rid_r0230 = 230
      integer, parameter :: rid_r0231 = 231
      integer, parameter :: rid_r0233 = 233
      integer, parameter :: rid_r0234 = 234
      integer, parameter :: rid_r0235 = 235
      integer, parameter :: rid_r0236 = 236
      integer, parameter :: rid_r0237 = 237
      integer, parameter :: rid_r0238 = 238
      integer, parameter :: rid_r0241 = 241
      integer, parameter :: rid_r0263 = 263
      integer, parameter :: rid_r0304 = 304
      integer, parameter :: rid_r0305 = 305
      end module M_RXT_ID_MOD
      module M_HET_ID_MOD
      implicit none
      end module M_HET_ID_MOD
