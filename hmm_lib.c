#include "util_lib.h"

void viterbi(HMM *hmm_ptr, char *O, char *output_buffer, char *aa_buffer,
             char *dna_buffer, char *sequence_head, bool whole_genome, bool format,
             int len_seq, char *dna, char *dna1, char *dna_f, char *dna_f1,
             char *protein, int *insert, int *c_delete, char *temp_str_ptr) {

    int *vpath;                          // optimal path after backtracking
    int **path;                          // viterbi path array
    double **alpha;                      // viterbi prob array
    int i, j, t, kk;
    int from, from0, to;   /*from0: i-2 position, from: i-1 position */
    int from2;             /* from2: i-2, i-1 for condition in probability */
    int gene_len;
    int num_d;          		/* the number of delete */
    double h_kd, r_kd, p_kd;
    double temp_alpha, prob;
    double start_freq;
    double final_score;

    int codon_start = 0;
    int dna_id = 0;
    int dna_f_id = 0;
    int out_nt;
    int start_t = -1;
    int end_t;
    int prev_match;
    int start_orf;
    int frame;
    int insert_id, delete_id;
    int temp_i[6]   = {0,0,0,0,0,0};
    int temp_i_1[6] = {1,1,1,1,1,1};
    int num_N = 0;


    /* Parse the sequence */
    Nucleotide sequence[len_seq];
    for (i = 0; i < len_seq; i++)
        sequence[i] = nt2int(O[i]);

    /***************************************************************/
    /* initialize                                                  */
    /***************************************************************/

    gene_len = (whole_genome)? 120 : 60;

    alpha = (double **)dmatrix(len_seq);
    path = (int **)imatrix(len_seq);
    vpath = (int *)ivector(len_seq);

    for (i = 0; i < NUM_STATE; i++) {
        alpha[i][0] = -1 * (hmm_ptr->pi[i]);
    }

    /* stop state */
    if ((sequence[0] == NUCL_T) && (((sequence[1] == NUCL_A) && (sequence[2] == NUCL_A)) ||
                                    ((sequence[1] == NUCL_A) && (sequence[2] == NUCL_G)) ||
                                    ((sequence[1] == NUCL_G) && (sequence[2] == NUCL_A)))) {

        alpha[E_STATE][0] = DBL_MAX;
        alpha[E_STATE][1] = DBL_MAX;
        path[E_STATE][1] = E_STATE;
        path[E_STATE][2] = E_STATE;

        alpha[M6_STATE][2] = DBL_MAX;
        alpha[M5_STATE][1] = DBL_MAX;
        alpha[M4_STATE][0] = DBL_MAX;
        alpha[M3_STATE][2] = DBL_MAX;
        alpha[M2_STATE][1] = DBL_MAX;
        alpha[M1_STATE][0] = DBL_MAX;

        if ((sequence[1] == NUCL_A) && (sequence[2] == NUCL_A)) {
            alpha[E_STATE][2] = alpha[E_STATE][2] - LOG_53;
        } else if ((sequence[1] == NUCL_A) && (sequence[2] == NUCL_G)) {
            alpha[E_STATE][2] = alpha[E_STATE][2] - LOG_16;
        } else if ((sequence[1] == NUCL_G) && (sequence[2] == NUCL_A)) {
            alpha[E_STATE][2] = alpha[E_STATE][2] - LOG_30;
        }
    }

    if ((sequence[2] == NUCL_A) &&
            (((sequence[0] == NUCL_T) && (sequence[1] == NUCL_T)) ||
             ((sequence[0] == NUCL_C) && (sequence[1] == NUCL_T)) ||
             ((sequence[0] == NUCL_T) && (sequence[1] == NUCL_C)))) {
        alpha[S_STATE_1][0] = DBL_MAX;
        alpha[S_STATE_1][1] = DBL_MAX;
        alpha[S_STATE_1][2] = alpha[S_STATE][0];
        path[S_STATE_1][1] = S_STATE_1;
        path[S_STATE_1][2] = S_STATE_1;

        alpha[M3_STATE_1][2] = DBL_MAX;
        alpha[M6_STATE_1][2] = DBL_MAX;

        if ((sequence[0] == NUCL_T) && (sequence[1] == NUCL_T)) {
            alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - LOG_53;
        } else if ((sequence[0] == NUCL_C) && (sequence[1] == NUCL_T)) {
            alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - LOG_16;
        } else if ((sequence[0] == NUCL_T) && (sequence[1] == NUCL_C)) {
            alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - LOG_30;
        }
    }

    int multiple=0;
    /******************************************************************/
    /*  fill out the rest of the columns                              */
    /******************************************************************/
    for (t = 1; t < len_seq; t++) {
        from = sequence[t-1];
        from0 = (t > 1)? sequence[t-2] : NUCL_G;
        to = sequence[t];

        /* if DNA is other than ACGT, do it later */
        if (from == 4) {
            from = 2;
        }
        if (from0 == 4) {
            from0 = 2;
        }
        if (to == 4) {
            to = 2;
            num_N += 1;
        } else {
            num_N = 0;
        }
        from2 = from0 * 4 + from;

        /******************/
        /* M state        */
        /******************/

        for (i = M1_STATE; i <= M6_STATE; i++)   {
            if (alpha[i][t]<DBL_MAX) {
                if (t==0) {
                } else {
                    if (i==M1_STATE) {
                        /* from M state */
                        j = M6_STATE;
                        alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_GG]) - (hmm_ptr->tr[TR_MM]) - (hmm_ptr->e_M[0][from2][to]);
                        path[i][t] = j;

                        /* from D state */
                        if (!whole_genome) {
                            for (j = M5_STATE; j >= M1_STATE; j--) {
                                if (j >= i ) {
                                    num_d = i-j+6;
                                } else if (j+1<i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if (num_d>0) {
                                    temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_MD]) -
                                                 (hmm_ptr->e_M[0][from2][to])
                                                 - LOG_25*(num_d-1) - (hmm_ptr->tr[TR_DD])*(num_d-2) -
                                                 (hmm_ptr->tr[TR_DM]);
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }

                        /* from Start state */
                        temp_alpha = alpha[S_STATE][t-1] - (hmm_ptr->e_M[0][from2][to]);
                        if ( temp_alpha < alpha[i][t] ) {
                            alpha[i][t] = temp_alpha;
                            path[i][t] = S_STATE;
                        }

                    } else {  /*i ==M2-M6*/

                        /* from M state */
                        j = i - 1;
                        alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_MM]) -
                                      (hmm_ptr->e_M[i-M1_STATE][from2][to]);
                        path[i][t] = j;


                        /* from D state */
                        if (!whole_genome) {
                            for (j = M6_STATE; j >= M1_STATE; j--) {
                                if (j >= i) {
                                    num_d = i-j+6;
                                } else if (j+1 < i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if (num_d > 0) {


                                    temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_MD]) -
                                                 (hmm_ptr->e_M[i-M1_STATE][from2][to])
                                                 - LOG_25*(num_d-1) - (hmm_ptr->tr[TR_DD])*(num_d-2) -
                                                 (hmm_ptr->tr[TR_DM]);
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }
                    }

                    /* from I state */
                    if (i == M1_STATE) {
                        j = I6_STATE;
                    } else {
                        j = I1_STATE + (i - M1_STATE -1);
                    }


                    /* to aviod stop codon */
                    if (t<2) {
                    } else if ((i==M2_STATE || i==M5_STATE) && (sequence[temp_i[j-I1_STATE]] == NUCL_T) &&
                               (((sequence[t] == NUCL_A) && (sequence[t+1] ==NUCL_A)) ||
                                ((sequence[t] == NUCL_A) && (sequence[t+1] ==NUCL_G)) ||
                                ((sequence[t] == NUCL_G) && (sequence[t+1] ==NUCL_A)))) {

                    } else if (((j-I1_STATE > 0) && (temp_i[j-I1_STATE]>0)) && ((i==M3_STATE || i==M6_STATE) && (sequence[temp_i[j-I1_STATE]-1] == NUCL_T) &&
                               (((sequence[temp_i[j-I1_STATE]] == NUCL_A) && (sequence[t] ==NUCL_A)) ||
                                ((sequence[temp_i[j-I1_STATE]] == NUCL_A) && (sequence[t] ==NUCL_G)) ||
                                ((sequence[temp_i[j-I1_STATE]] == NUCL_G) && (sequence[t] ==NUCL_A))))) {
                    } else {
                        temp_alpha = alpha[j][t-1]  - (hmm_ptr->tr[TR_IM]) - LOG_25;
                        if ( temp_alpha < alpha[i][t]) {
                            alpha[i][t] = temp_alpha;
                            path[i][t] = j;
                        }
                    }
                }
            }
        }

        /******************/
        /* I state        */
        /******************/
        for (i = I1_STATE; i <= I6_STATE; i++) {

            if (t != 0) {
                /* from I state */
                j = i;
                alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_II]) -
                              (hmm_ptr->tr_I_I[from][to]);
                path[i][t] = j;

                /* from M state */
                j = i - I1_STATE + M1_STATE ;
                if (i == I6_STATE) {
                    temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_GG]) -
                                 (hmm_ptr->tr[TR_MI]) -(hmm_ptr->tr_M_I[from][to]);
                } else {
                    temp_alpha = alpha[j][t-1]  -
                                 (hmm_ptr->tr[TR_MI]) -(hmm_ptr->tr_M_I[from][to]);
                }
                if (temp_alpha < alpha[i][t]) {
                    alpha[i][t] = temp_alpha;
                    path[i][t] = j;

                    temp_i[i-I1_STATE] = t-1;
                }
            }
        }

        /******************/
        /* M' state        */
        /******************/

        for (i = M1_STATE_1; i <= M6_STATE_1; i++)   {
            if  ((i==M1_STATE_1 || i==M4_STATE_1)&& t>=3 &&
                    (((sequence[t-3] == NUCL_T) && (sequence[t-2] == NUCL_T) && (sequence[t-1] == NUCL_A)) ||
                     ((sequence[t-3] == NUCL_C) && (sequence[t-2] == NUCL_T) && (sequence[t-1] == NUCL_A)) ||
                     ((sequence[t-3] == NUCL_T) && (sequence[t-2] == NUCL_C) && (sequence[t-1] == NUCL_A)))) {

                /* from Start state  since this is actually stop codon in minus strand */
                alpha[i][t] = alpha[S_STATE_1][t-1] -
                              (hmm_ptr->e_M_1[i-M1_STATE_1][from2][to]);
                path[i][t] = S_STATE_1;

            } else {
                if (t != 0) {
                    if (i == M1_STATE_1 ) {

                        /* from M state */
                        j = M6_STATE_1;
                        alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_GG]) -
                                      (hmm_ptr->tr[TR_MM]) - (hmm_ptr->e_M_1[0][from2][to]);
                        path[i][t] = j;

                        /* from D state */
                        if (!whole_genome) {
                            for (j = M5_STATE_1; j >= M1_STATE_1; j--) {
                                if (j >= i) {
                                    num_d = i-j+6;
                                } else if (j+1 <i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if (num_d > 0) {
                                    temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_MD]) -
                                                 (hmm_ptr->e_M_1[0][from2][to])
                                                 - LOG_25*(num_d-1) - (hmm_ptr->tr[TR_DD])*(num_d-2) -
                                                 (hmm_ptr->tr[TR_DM]);
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }

                    } else {

                        /* from M state */
                        j = i - 1;
                        alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_MM]) -
                                      (hmm_ptr->e_M_1[i-M1_STATE_1][from2][to]);
                        path[i][t] = j;

                        /* from D state */
                        if (!whole_genome) {
                            for (j = M6_STATE_1; j>=M1_STATE_1; j--) {
                                if (j >= i ) {
                                    num_d = i-j+6;
                                } else if (j+1 < i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if (num_d>0) {
                                    temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_MD]) -
                                                 (hmm_ptr->e_M_1[i-M1_STATE_1][from2][to])
                                                 - LOG_25*(num_d-1) - (hmm_ptr->tr[TR_DD])*(num_d-2) -
                                                 (hmm_ptr->tr[TR_DM]);
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }
                    }

                    /* from I state */
                    if (i == M1_STATE_1) {
                        j = I6_STATE_1;
                    } else {
                        j = I1_STATE_1 + (i - M1_STATE_1 -1);
                    }


                    //!! 11 --> 16
                    //!! What is the actual point of -1

                    /* to aviod stop codon */
                    if (t<2) {
                    } else  if ((i==M2_STATE_1 ||
                                 i==M5_STATE_1) &&
                                (sequence[t+1] == NUCL_A ) &&
                                (((sequence[temp_i_1[j-I1_STATE_1]] == NUCL_T) && (sequence[t] ==NUCL_T)) ||
                                 ((sequence[temp_i_1[j-I1_STATE_1]] == NUCL_C) && (sequence[t] ==NUCL_T)) ||
                                 ((sequence[temp_i_1[j-I1_STATE_1]] == NUCL_T) && (sequence[t] ==NUCL_C)))) {

                    } else if ((i==M3_STATE_1 ||
                                i==M6_STATE_1) &&
                               (sequence[t] == NUCL_A ) &&
                               (((sequence[temp_i_1[j-I1_STATE_1]-1] == NUCL_T) && (sequence[temp_i_1[j-I1_STATE_1]] ==NUCL_T)) ||
                                ((sequence[temp_i_1[j-I1_STATE_1]-1] == NUCL_C) && (sequence[temp_i_1[j-I1_STATE_1]] ==NUCL_T)) ||
                                ((sequence[temp_i_1[j-I1_STATE_1]-1] == NUCL_T) && (sequence[temp_i_1[j-I1_STATE_1]] ==NUCL_C)))) {
                    } else {

                        temp_alpha = alpha[j][t-1]  - (hmm_ptr->tr[TR_IM]) - LOG_25;
                        if ( temp_alpha < alpha[i][t]) {
                            alpha[i][t] = temp_alpha;
                            path[i][t] = j;
                        }
                    }
                }
            }
        }

        /******************/
        /* I' state        */
        /******************/
        for (i = I1_STATE_1; i <= I6_STATE_1; i++) {
            if (t != 0) {
                /* from I state */
                j = i;
                alpha[i][t] = alpha[j][t-1] - (hmm_ptr->tr[TR_II]) -
                              (hmm_ptr->tr_I_I[from][to]);
                path[i][t] = j;

                if (t<5) continue;
                /* from M state */
                if (path[S_STATE_1][t-3]!= R_STATE && path[S_STATE_1][t-4] !=R_STATE &&
                        path[S_STATE_1][t-5] !=R_STATE) {
                    j = i - I1_STATE_1 + M1_STATE_1;
                    if (i==I6_STATE_1) {
                        temp_alpha = alpha[j][t-1] - (hmm_ptr->tr[TR_GG]) -
                                     (hmm_ptr->tr[TR_MI]) -(hmm_ptr->tr_M_I[from][to]);
                    } else {
                        temp_alpha = alpha[j][t-1]  -
                                     (hmm_ptr->tr[TR_MI]) -(hmm_ptr->tr_M_I[from][to]);
                    }
                    if (temp_alpha < alpha[i][t]) {
                        alpha[i][t] = temp_alpha;
                        path[i][t] = j;
                        //!! We are addressing the character array with this.
                        temp_i_1[i-I1_STATE_1] = t-1;
                    }
                }
            }
        }

        /***********************/
        /* Non_coding state    */
        /***********************/

        if (t != 0) {
            alpha[R_STATE][t] = alpha[R_STATE][t-1] - (hmm_ptr->tr_R_R[from][to]) -  (hmm_ptr->tr[TR_RR]);
            path[R_STATE][t] = R_STATE;

            temp_alpha = alpha[E_STATE][t-1]  - (hmm_ptr->tr[TR_ER])  ;
            if (temp_alpha < alpha[R_STATE][t] ) {
                alpha[R_STATE][t] = temp_alpha;
                path[R_STATE][t] = E_STATE;
            }

            temp_alpha = alpha[E_STATE_1][t-1] - (hmm_ptr->tr[TR_ER]) ;
            if (temp_alpha < alpha[R_STATE][t] ) {
                alpha[R_STATE][t] = temp_alpha;
                path[R_STATE][t] = E_STATE_1;
            }
            alpha[R_STATE][t] -= LOG_95;
        }

        /******************/
        /* END state      */
        /******************/
        if (alpha[E_STATE][t] == 0) {

            alpha[E_STATE][t] = DBL_MAX;
            path[E_STATE][t] = NOSTATE;

            if (t < len_seq -2 && (sequence[t] == NUCL_T)  &&
                    (((sequence[t+1] == NUCL_A) && (sequence[t+2] == NUCL_A)) ||
                     ((sequence[t+1] == NUCL_A) && (sequence[t+2] == NUCL_G)) ||
                     ((sequence[t+1] == NUCL_G) && (sequence[t+2] == NUCL_A)))) {

                alpha[E_STATE][t+2] = DBL_MAX;
                /* transition from frame4,frame5,and frame6 */
                temp_alpha = alpha[M6_STATE][t-1] - (hmm_ptr->tr[TR_GE]);
                if (temp_alpha < alpha[E_STATE][t+2]) {
                    alpha[E_STATE][t+2] = temp_alpha;
                    path[E_STATE][t] = M6_STATE;
                }

                /* transition from frame1,frame2,and frame3 */
                temp_alpha  = alpha[M3_STATE][t-1] - (hmm_ptr->tr[TR_GE]);
                if (temp_alpha < alpha[E_STATE][t+2]) {
                    alpha[E_STATE][t+2] = temp_alpha;
                    path[E_STATE][t] = M3_STATE;
                }

                alpha[E_STATE][t] = DBL_MAX;
                alpha[E_STATE][t+1] = DBL_MAX;
                path[E_STATE][t+1] = E_STATE;
                path[E_STATE][t+2] = E_STATE;

                alpha[M6_STATE][t+2] = DBL_MAX;
                alpha[M5_STATE][t+1] = DBL_MAX;
                alpha[M4_STATE][t] = DBL_MAX;
                alpha[M3_STATE][t+2] = DBL_MAX;
                alpha[M2_STATE][t+1] = DBL_MAX;
                alpha[M1_STATE][t] = DBL_MAX;

                if ((sequence[t+1] == NUCL_A) && (sequence[t+2] == NUCL_A)) {
                    alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - (0.54);
                } else if ((sequence[t+1] == NUCL_A) && (sequence[t+2] == NUCL_G)) {
                    alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - LOG_16;
                } else if ((sequence[t+1] == NUCL_G) && (sequence[t+2] == NUCL_A)) {
                    alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - LOG_30;
                }

                /* adjustment based on probability distribution */
                start_freq=0;

                double sub_sum = 0;

                if (t>=60) { /* bug reported by Yu-Wei */
                    for (i = -60; i <= -3; i++) {
                        if (t+i+2 < len_seq) {
                            start_freq -= (hmm_ptr->tr_E[i+60][trinucleotide(sequence[t+i], sequence[t+i+1], sequence[t+i+2])]);
                        }
                    }
                } else {
                    for (i = (-1*t); i <= -3; i++) {
                        if (t+i+2 < len_seq) {
                            sub_sum += (hmm_ptr->tr_E[i+60][trinucleotide(sequence[t+i], sequence[t+i+1], sequence[t+i+2])]);
                        }
                    }
                    sub_sum = sub_sum * 58 / (-3 + t + 1);
                    start_freq -= sub_sum;
                }

                h_kd = hmm_ptr->E_dist[2] * exp(-1*pow(start_freq-hmm_ptr->E_dist[1],2)/(2*pow(hmm_ptr->E_dist[0],2)));
                r_kd = hmm_ptr->E_dist[5] * exp(-1*pow(start_freq-hmm_ptr->E_dist[4],2)/(2*pow(hmm_ptr->E_dist[3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log(p_kd);
            }
        }

        /*************************************************/
        /* START' state                                  */
        /* origianlly stop codon of genes in - strand    */
        /*************************************************/
        if (alpha[S_STATE_1][t] == 0) {

            alpha[S_STATE_1][t] = DBL_MAX;
            path[S_STATE_1][t] = NOSTATE;


            if (t<len_seq-2 && (sequence[t+2] == NUCL_A) &&
                    (((sequence[t] == NUCL_T) && (sequence[t+1] == NUCL_T)) ||
                     ((sequence[t] == NUCL_C) && (sequence[t+1] == NUCL_T)) ||
                     ((sequence[t] == NUCL_T) && (sequence[t+1] == NUCL_C)))) {

                alpha[S_STATE_1][t] = DBL_MAX;
                path[S_STATE_1][t] = R_STATE;
                alpha[S_STATE_1][t+1] = DBL_MAX;
                alpha[S_STATE_1][t+2] = alpha[R_STATE][t-1] - (hmm_ptr->tr[TR_RS]);
                path[S_STATE_1][t+1] = S_STATE_1;
                path[S_STATE_1][t+2] = S_STATE_1;

                temp_alpha = alpha[E_STATE_1][t-1] - (hmm_ptr->tr[TR_ES]);
                if (temp_alpha < alpha[S_STATE_1][t+2]) {
                    alpha[S_STATE_1][t+2] = temp_alpha;
                    path[S_STATE_1][t] = E_STATE_1;
                }

                temp_alpha = alpha[E_STATE][t-1] - (hmm_ptr->tr[TR_ES1]);
                if (temp_alpha < alpha[S_STATE_1][t+2]) {
                    alpha[S_STATE_1][t+2] = temp_alpha;
                    path[S_STATE_1][t] = E_STATE;
                }

                alpha[M3_STATE_1][t+2] = DBL_MAX;
                alpha[M6_STATE_1][t+2] = DBL_MAX;

                if ((sequence[t] == NUCL_T) && (sequence[t+1] == NUCL_T)) {
                    alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - (0.54);
                } else if ((sequence[t] == NUCL_C) && (sequence[t+1] == NUCL_T)) {
                    alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - LOG_16;
                } else if ((sequence[t] == NUCL_T) && (sequence[t+1] == NUCL_C)) {
                    alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - LOG_30;
                }

                /* adjustment based on probability distribution */
                start_freq=0;
                for (i = 3; i <= 60; i++) {
                    if (t+i+2 < len_seq) {
                        start_freq -= (hmm_ptr->tr_S_1[i-3][trinucleotide(sequence[t+i], sequence[t+i+1], sequence[t+i+2])]);
                    }
                }
                h_kd = hmm_ptr->S1_dist[2] * exp(-1*pow(start_freq-hmm_ptr->S1_dist[1],2)/(2*pow(hmm_ptr->S1_dist[0],2)));
                r_kd = hmm_ptr->S1_dist[5] * exp(-1*pow(start_freq-hmm_ptr->S1_dist[4],2)/(2*pow(hmm_ptr->S1_dist[3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log(p_kd);
            }
        }

        /************************/
        /* START state          */
        /************************/
        if (alpha[S_STATE][t] == 0) {

            alpha[S_STATE][t] = DBL_MAX;
            path[S_STATE][t] = NOSTATE;

            if (t<len_seq-2 &&  (sequence[t+1] == NUCL_T) && (sequence[t+2] == NUCL_G)&&
                    ((sequence[t] == NUCL_A) || (sequence[t] == NUCL_G) ||  (sequence[t] == NUCL_T))) {

                alpha[S_STATE][t] = DBL_MAX;
                alpha[S_STATE][t+1] = DBL_MAX;
                alpha[S_STATE][t+2] = alpha[R_STATE][t-1] - (hmm_ptr->tr[TR_RS]);
                path[S_STATE][t] = R_STATE;
                path[S_STATE][t+1] = S_STATE;
                path[S_STATE][t+2] = S_STATE;

                temp_alpha = alpha[E_STATE][t-1] - (hmm_ptr->tr[TR_ES]);
                if (temp_alpha < alpha[S_STATE][t+2]) {
                    alpha[S_STATE][t+2] = temp_alpha;
                    path[S_STATE][t] = E_STATE;
                }

                temp_alpha = alpha[E_STATE_1][t-1] - (hmm_ptr->tr[TR_ES1]);
                if (temp_alpha < alpha[S_STATE][t+2]) {
                    alpha[S_STATE][t+2] = temp_alpha;
                    path[S_STATE][t] = E_STATE_1;
                }


                if (sequence[t] == NUCL_A) {
                    alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - LOG_83;
                } else if (sequence[t] == NUCL_G) {
                    alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - LOG_10;
                } else if (sequence[t] == NUCL_T) {
                    alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - LOG_07;
                }

                /* adjustment based on probability distribution */
                start_freq=0;

                double sub_sum = 0;

                if (t>=30) {
                    for (i = -30; i <= 30; i++) {
                        if (t+i+2 < len_seq) {
                            start_freq -= (hmm_ptr->tr_S[i+30][trinucleotide(sequence[t+i], sequence[t+i+1], sequence[t+i+2])]);
                        }
                    }
                } else {
                    for (i = (-1*t); i <= 30; i++) {
                        if (t+i+2 < len_seq) {
                            sub_sum += (hmm_ptr->tr_S[i+30][trinucleotide(sequence[t+i], sequence[t+i+1], sequence[t+i+2])]);
                        }
                    }
                    sub_sum = sub_sum * 61 / (30 + t + 1);
                    start_freq -= sub_sum;
                }

                h_kd = hmm_ptr->S_dist[2] * exp(-1*pow(start_freq-hmm_ptr->S_dist[1],2)/(2*pow(hmm_ptr->S_dist[0],2)));
                r_kd = hmm_ptr->S_dist[5] * exp(-1*pow(start_freq-hmm_ptr->S_dist[4],2)/(2*pow(hmm_ptr->S_dist[3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log(p_kd);

            }
        }

        /**********************************************/
        /* END' state                                 */
        /* origianlly start codon of genes in - strand */
        /**********************************************/
        if (alpha[E_STATE_1][t] == 0) {

            alpha[E_STATE_1][t] = DBL_MAX;
            path[E_STATE_1][t] = NOSTATE;

            if (t < len_seq - 2 && (sequence[t] == NUCL_C) && (sequence[t+1] == NUCL_A) &&
                    ((sequence[t+2] == NUCL_T) || (sequence[t+2] == NUCL_C) || (sequence[t+2] == NUCL_A))) {

                /* transition from frame6 */
                alpha[E_STATE_1][t+2] = alpha[M6_STATE_1][t-1] - (hmm_ptr->tr[TR_GE]);
                path[E_STATE_1][t] = M6_STATE_1;
                alpha[E_STATE_1][t] = DBL_MAX;
                alpha[E_STATE_1][t+1] = DBL_MAX;
                path[E_STATE_1][t+1] = E_STATE_1;
                path[E_STATE_1][t+2] = E_STATE_1;

                if (sequence[t+2] == NUCL_T) {
                    alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - LOG_83;
                } else if (sequence[t+2] == NUCL_C ) {
                    alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - LOG_10;
                } else if (sequence[t+2] == NUCL_A ) {
                    alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - LOG_07;
                }

                /* adjustment based on probability distribution */
                start_freq=0;

                double sub_sum = 0;

                if (t>=30) {
                    for (i = -30; i <= 30; i++) {
                        if (t+i+2 < len_seq) {
                            start_freq -= (hmm_ptr->tr_E_1[i+30][trinucleotide(sequence[t+i], sequence[t+i+1], sequence[t+i+2])]);
                        }
                    }
                } else {
                    for (i = (-1*t); i <= 30; i++) {
                        if (t+i+2 < len_seq) {
                            sub_sum += (hmm_ptr->tr_E_1[i+30][trinucleotide(sequence[t+i], sequence[t+i+1], sequence[t+i+2])]);
                        }
                    }
                    sub_sum = sub_sum * 61 / (30 + t + 1);
                    start_freq -= sub_sum;
                }

                h_kd = hmm_ptr->E1_dist[2] *
                       exp(-1*pow(start_freq-hmm_ptr->E1_dist[1],2)/(2*pow(hmm_ptr->E1_dist[0],2)))
                       ;
                r_kd = hmm_ptr->E1_dist[5] *
                       exp(-1*pow(start_freq-hmm_ptr->E1_dist[4],2)/(2*pow(hmm_ptr->E1_dist[3],2)))
                       ;
                p_kd = h_kd / (h_kd + r_kd);

                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log(p_kd);
            }
        }
        if (num_N>9) {

            for (i = 0; i < NUM_STATE; i++) {
                if (i!=R_STATE) {
                    alpha[i][t] = DBL_MAX;
                    path[i][t] = R_STATE;
                }
            }
        }
    }




    /***********************************************************/
    /* backtrack array to find the optimal path                */
    /***********************************************************/


    sprintf(output_buffer, "%s\n", sequence_head);

    /* find the state for sequence[N] with the highest probability */
    prob = DBL_MAX;
    for (i = 0; i < NUM_STATE; i++) {

        if (alpha[i][len_seq-1] < prob) {
            prob = alpha[i][len_seq-1];
            vpath[len_seq-1] = i;
        }
    }

    /* backtrack the opitmal path */
    for (t = len_seq-2; t >= 0; t--) {
        if (t+1 < 0 || vpath[t+1] < 0) continue;
        vpath[t] = path[vpath[t+1]][t+1];
    }

    for (t = 0; t < len_seq; t++) {

        if (codon_start == 0 && start_t < 0 &&
                ((vpath[t]>=M1_STATE && vpath[t]<=M6_STATE) ||
                 (vpath[t]>=M1_STATE_1 && vpath[t]<=M6_STATE_1) ||
                 vpath[t] == S_STATE || vpath[t] == S_STATE_1 )) {
            start_t=t+1;
        }

        if (codon_start == 0 &&
                (vpath[t]==M1_STATE || vpath[t]==M4_STATE ||
                 vpath[t]==M1_STATE_1 || vpath[t]==M4_STATE_1)) {

            stopMemset(dna, STRINGLEN);
            stopMemset(dna1, STRINGLEN);//
            stopMemset(dna_f, STRINGLEN);//
            stopMemset(dna_f1, STRINGLEN);//
            stopMemset(protein, STRINGLEN);
            stopMemset(insert, STRINGLEN);//
            stopMemset(c_delete, STRINGLEN);//

            insert_id = 0;
            delete_id = 0;
            dna_id = 0;
            dna_f_id = 0;
            dna[dna_id] = O[t];
            dna_f[dna_f_id] = O[t];
            start_orf = t+1;
            prev_match = vpath[t];

            if (vpath[t] < M6_STATE) {
                codon_start = 1;
            } else {
                codon_start = -1;
            }

        } else if (codon_start != 0 && (vpath[t] == E_STATE || vpath[t] == E_STATE_1 || t == len_seq-1)) {

            if (vpath[t] == E_STATE || vpath[t] == E_STATE_1) {
                end_t = t+3;
            } else {
                end_t = t+1;

                /* FGS1.12 start: remove incomplete codon */
                int temp_t = t;
                while (vpath[temp_t] != M1_STATE && vpath[temp_t] != M4_STATE  &&
                        vpath[temp_t] != M1_STATE_1  && vpath[temp_t] != M4_STATE_1) {
                    dna_f[dna_f_id] = '\0';
                    dna_f_id--;

                    dna[dna_id] = '\0';
                    dna_id--;

                    temp_t--;
                }
                /* FGS1.12 end: remove incomplete codon */
            }
            final_score = (alpha[vpath[end_t-4]][end_t-4] - alpha[vpath[start_t+2]][start_t+2] )/(end_t-start_t-5);
            frame = start_orf%3;
            if (frame==0) {
                frame=3;
            }

            //!! Transfer all of the output buffer writing code to another function. Modularize this.

            if (dna_id > gene_len) {
                print_outputs(codon_start, start_t, end_t, frame, output_buffer, aa_buffer, dna_buffer, sequence_head,
                              dna, dna_id + 1, dna1, dna_f, dna_f1, protein, insert, c_delete, insert_id, delete_id, format, temp_str_ptr,multiple);
                multiple++;
            }

            codon_start = 0;
            start_t = -1;
            end_t = -1;
            dna_id = 0;
            dna_f_id = 0;

        } else if (codon_start != 0 &&
                   ((vpath[t] >= M1_STATE && vpath[t] <= M6_STATE) ||
                    (vpath[t] >= M1_STATE_1 && vpath[t] <= M6_STATE_1)) &&
                   vpath[t] - prev_match < 6) {

            if (vpath[t] < prev_match) {
                out_nt = vpath[t]+6-prev_match;
            } else {
                out_nt = vpath[t]-prev_match;
            }
            for (kk = 0; kk < out_nt; kk++) {  /* for deleted nt in reads */
                dna_id ++;
                dna[dna_id] = 'N';
                dna_f_id++;
                dna_f[dna_f_id] = 'x';
                if (kk>0) {
                    c_delete[delete_id] = t+1;
                    delete_id++;
                }
            }
            dna[dna_id] = O[t];
            dna_f[dna_f_id] = O[t];
            prev_match = vpath[t];

        } else if (codon_start != 0 &&
                   ((vpath[t] >= I1_STATE && vpath[t] <= I6_STATE) ||
                    (vpath[t] >= I1_STATE_1 && vpath[t] <= I6_STATE_1))) {
            dna_f_id ++;
            dna_f[dna_f_id] = tolower(sequence[t]);
            insert[insert_id] = t+1;
            insert_id++;

        } else if (codon_start != 0 && vpath[t] == R_STATE) {
            /* for long NNNNNNNNN, pretend R state */
            codon_start = 0;
            start_t = -1;
            end_t = -1;
            dna_id = 0;
            dna_f_id = 0;

        }
    }
    free_dmatrix(alpha);
    free_imatrix(path);
    free(vpath);

    vpath = 0;
    dna = 0;
    dna1 = 0;
    dna_f = 0;
    dna_f = 0;
    protein = 0;
}

void get_prob_from_cg(HMM *hmm_ptr, TRAIN *train_ptr, char *O, int len_seq) {

    int cg_id = -1;
    int cg_count=0;
    int i,j,k;

    for (i = 0; i < len_seq; i++) {
        if ((O[i] == 'C'||O[i] =='c') || (O[i] == 'G'||O[i] == 'g') ) {
            cg_count++;
        }
    }

    cg_count = floor((cg_count*1.0/len_seq)*100)-26;

    if (cg_count < 0) {
        cg_count = 0;
    } else if (cg_count > 43) {
        cg_count = 43;
    }

    memcpy(hmm_ptr->e_M, train_ptr->trans[cg_count], sizeof(hmm_ptr->e_M));
    memcpy(hmm_ptr->e_M_1, train_ptr->rtrans[cg_count], sizeof(hmm_ptr->e_M_1));
    memcpy(hmm_ptr->tr_R_R, train_ptr->noncoding[cg_count], sizeof(hmm_ptr->tr_R_R));
    memcpy(hmm_ptr->tr_S, train_ptr->start[cg_count], sizeof(hmm_ptr->tr_S));
    memcpy(hmm_ptr->tr_E, train_ptr->stop[cg_count], sizeof(hmm_ptr->tr_E));
    memcpy(hmm_ptr->tr_S_1, train_ptr->start1[cg_count], sizeof(hmm_ptr->tr_S_1));
    memcpy(hmm_ptr->tr_E_1, train_ptr->stop1[cg_count], sizeof(hmm_ptr->tr_E_1));
    memcpy(hmm_ptr->S_dist, train_ptr->S_dist[cg_count], sizeof(hmm_ptr->S_dist));
    memcpy(hmm_ptr->E_dist, train_ptr->E_dist[cg_count], sizeof(hmm_ptr->E_dist));
    memcpy(hmm_ptr->S1_dist, train_ptr->S1_dist[cg_count], sizeof(hmm_ptr->S1_dist));
    memcpy(hmm_ptr->E1_dist, train_ptr->E1_dist[cg_count], sizeof(hmm_ptr->E1_dist));
}


//!! This function is for some strange reason definetly leaking memory. Need to find out why.
void get_train_from_file(char *filename, HMM *hmm_ptr, char *mfilename,
                         char *mfilename1, char *nfilename, char *sfilename, char *pfilename,
                         char *s1filename, char *p1filename, char *dfilename, TRAIN *train_ptr) {

    int i, j, k, p;
    double prob;
    FILE *fp, *fpm, *fpm1, *fpn, *fps, *fpp, *fps1, *fpp1, *fpd;

    char name[10];
    char head[20];
    char start[10];
    char end[10];

    /****************************************************/
    /* transition                                       */
    /****************************************************/
    fp = fopen (filename , "r");

    /* Transition */
    fscanf(fp, "%s", head);
    for (i = 0; i < NUM_TRANSITIONS; i++) {
        //!! This causes a memory leak for unknown reasons.
        fscanf(fp, "%s %lf", name, &prob);
        hmm_ptr->tr[tr2int(name)] = log(prob);
    }

    /* TransitionMI */
    fscanf(fp, "%s", head);
    for (i = 0; i < 16; i++) {
        fscanf(fp, "%s %s %lf\n", start, end, &prob);
        hmm_ptr->tr_M_I[nt2int(start[0])][nt2int(end[0])] = log(prob);
    }

    /* TransitionII */
    fscanf(fp, "%s", head);
    for (i = 0; i < 16; i++) {
        fscanf(fp, "%s %s %lf", start, end, &prob);
        hmm_ptr->tr_I_I[nt2int(start[0])][nt2int(end[0])] = log(prob);
    }

    /* PI */
    fscanf(fp, "%s", head);
    for (i = 0; i < NUM_STATE; i++) {
        fscanf(fp, "%s %lf", name, &prob);
        hmm_ptr->pi[i] = log(prob);
    }
    fclose(fp);

    /****************************************************/
    /* M state transition                               */
    /****************************************************/
    fpm = fopen (mfilename , "r");
    for (p = 0; p < 44; p++) {                       /* cg */
        fscanf(fpm, "%s", head);
        for (i = 0; i < 6; i++) {                      /* period */
            for (j = 0; j < 16; j++) {                   /* condition */
                for (k = 0; k < 4; k++) {                  /* emission */
                    fscanf(fpm, "%lf", &prob);
                    train_ptr->trans[p][i][j][k] = log(prob);
                }
            }
        }
    }
    fclose(fpm);

    /****************************************************/
    /* M state_1 transition                             */
    /****************************************************/
    fpm1 = fopen (mfilename1 , "r");
    for (p = 0; p < 44; p++) {
        fscanf(fpm1, "%s", head);
        for (i = 0; i < 6; i++) {
            for (j = 0; j < 16; j++) {
                for (k = 0; k < 4; k++) {
                    fscanf(fpm1, "%lf", &prob);
                    train_ptr->rtrans[p][i][j][k] = log(prob);
                }
            }
        }
    }
    fclose(fpm1);

    /****************************************************/
    /* noncoding state  transition                      */
    /****************************************************/
    fpn = fopen (nfilename, "r");
    for (p = 0; p < 44; p++) {
        fscanf(fpn, "%s", head);
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {
                fscanf(fpn, "%lf", &prob);
                train_ptr->noncoding[p][j][k] = log(prob);
            }
        }
    }
    fclose(fpn);

    /****************************************************/
    /* start                                            */
    /****************************************************/
    fps = fopen (sfilename, "r");
    for (p = 0; p < 44; p++) {
        fscanf(fps, "%s", head);
        for (j = 0; j < 61; j++) {
            for (k = 0; k < 64; k++) {
                fscanf(fps, "%lf", &prob);
                train_ptr->start[p][j][k] = log(prob);
            }
        }
    }
    fclose(fps);

    /****************************************************/
    /* stop                                             */
    /****************************************************/
    fpp = fopen (sfilename, "r");
    for (p = 0; p < 44; p++) {
        fscanf(fpp, "%s", head);
        for (j = 0; j < 58; j++) {
            for (k = 0; k < 64; k++) {
                fscanf(fpp, "%lf", &prob);
                train_ptr->stop[p][j][k] = log(prob);
            }
        }
    }
    fclose(fpp);

    /****************************************************/
    /* start1                                           */
    /****************************************************/
    fps1 = fopen (s1filename, "r");
    for (p = 0; p < 44; p++) {
        fscanf(fps1, "%s", head);
        for (j = 0; j < 58; j++) {
            for (k = 0; k < 64; k++) {
                fscanf(fps1, "%lf", &prob);
                train_ptr->start1[p][j][k] = log(prob);
            }
        }
    }
    fclose(fps1);

    /****************************************************/
    /* stop1                                            */
    /****************************************************/
    fpp1 = fopen (p1filename, "r");
    for (p = 0; p < 44; p++) {
        fscanf(fpp1, "%s", head);
        for (j = 0; j < 61; j++) {
            for (k = 0; k < 64; k++) {
                fscanf(fpp1, "%lf", &prob);
                train_ptr->stop1[p][j][k] = log(prob);
            }
        }
    }
    fclose(fpp1);


    /****************************************************/
    /* pwm distribution                                 */
    /****************************************************/
    fpd = fopen (dfilename, "r");
    for (p = 0; p < 44; p++) {
        fscanf(fpd, "%s", head);
        for (k = 0; k < 6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->S_dist[p][k] = prob;
        }
        for (k = 0; k < 6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->E_dist[p][k] = prob;
        }
        for (k = 0; k < 6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->S1_dist[p][k] = prob;
        }
        for (k = 0; k < 6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->E1_dist[p][k] = prob;
        }
    }
    fclose(fpd);

}

void print_outputs(int codon_start, int start_t, int end_t, int frame, char *output_buffer, char *aa_buffer, char *dna_buffer,
                   char *sequence_head_short, char *dna, int dna_len, char *rc_dna, char *dna_f, char *rc_dna_f, char *protein,
                   int *insertions, int *deletions, int insertions_len, int deletions_len, bool format, char *temp_str_ptr, unsigned int multiple) {
    int i;
    char strand_sign = (codon_start == 1)? '+' : '-';

    if (codon_start != 1 && codon_start != -1)
        return;

    /* Print the insertions and deletions to the output buffer */
    sprintf(temp_str_ptr, "%d\t%d\t%c\t%d\t", start_t, end_t, strand_sign, frame);
    strcat(output_buffer, temp_str_ptr);
    strcat(output_buffer, "I:");
    for (i = 0; i < insertions_len; i++) {
        sprintf(temp_str_ptr, "%d,", insertions[i]);
        strcat(output_buffer, temp_str_ptr);
    }
    strcat(output_buffer, "\tD:");
    for (i = 0; i < deletions_len; i++) {
        sprintf(temp_str_ptr, "%d,", deletions[i]);
        strcat(output_buffer, temp_str_ptr);
    }
    strcat(output_buffer, "\n");

    /* Now fill the AA buffer with the translated proteins */
    if (multiple)
        strcat(aa_buffer, "\t");
    sprintf(temp_str_ptr, "%s_%d_%d_%c\n", sequence_head_short, start_t, end_t, strand_sign);
    strcat(aa_buffer, temp_str_ptr);
    get_protein(dna, dna_len, protein, codon_start);
    sprintf(temp_str_ptr, "%s\n", protein);
    strcat(aa_buffer, temp_str_ptr);

    /* Similar for the DNA buffer */
    sprintf(temp_str_ptr, "%s_%d_%d_%c\n", sequence_head_short, start_t, end_t, strand_sign);
    strcat(dna_buffer, temp_str_ptr);
    /* Don't forget to print the reverse complement if in opposite strand */
    if (codon_start == 1) {
        sprintf(temp_str_ptr, "%s\n", (format)? dna_f : dna);
    } else {
        get_rc_dna(dna, dna_len, rc_dna);
        get_rc_dna_indel(dna_f, dna_len, rc_dna_f);
        sprintf(temp_str_ptr, "%s\n", (format)? rc_dna_f : rc_dna);
    }
    strcat(dna_buffer, temp_str_ptr);
}
