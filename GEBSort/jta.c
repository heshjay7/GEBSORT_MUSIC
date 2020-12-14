
/* collection of code related to jta */

#include "gdecomp.h"

#include "GEBSort.h"
#include "GTMerge.h"

//useful bit mask values (JTA)

const unsigned short usCLR_BIT0 = 0xFFFE;
const unsigned short usCLR_BIT1 = 0xFFFD;
const unsigned short usCLR_BIT2 = 0xFFFB;
const unsigned short usCLR_BIT3 = 0xFFF7;
const unsigned short usCLR_BIT4 = 0xFFEF;
const unsigned short usCLR_BIT5 = 0xFFDF;
const unsigned short usCLR_BIT6 = 0xFFBF;
const unsigned short usCLR_BIT7 = 0xFF7F;
const unsigned short usCLR_BIT8 = 0xFEFF;
const unsigned short usCLR_BIT9 = 0xFDFF;
const unsigned short usCLR_BIT10 = 0xFBFF;
const unsigned short usCLR_BIT11 = 0xF7FF;
const unsigned short usCLR_BIT12 = 0xEFFF;
const unsigned short usCLR_BIT13 = 0xDFFF;
const unsigned short usCLR_BIT14 = 0xBFFF;
const unsigned short usCLR_BIT15 = 0x7FFF;

const unsigned short usSET_BIT0 = 0x0001;
const unsigned short usSET_BIT1 = 0x0002;
const unsigned short usSET_BIT2 = 0x0004;
const unsigned short usSET_BIT3 = 0x0008;
const unsigned short usSET_BIT4 = 0x0010;
const unsigned short usSET_BIT5 = 0x0020;
const unsigned short usSET_BIT6 = 0x0040;
const unsigned short usSET_BIT7 = 0x0080;
const unsigned short usSET_BIT8 = 0x0100;
const unsigned short usSET_BIT9 = 0x0200;
const unsigned short usSET_BIT10 = 0x0400;
const unsigned short usSET_BIT11 = 0x0800;
const unsigned short usSET_BIT12 = 0x1000;
const unsigned short usSET_BIT13 = 0x2000;
const unsigned short usSET_BIT14 = 0x4000;
const unsigned short usSET_BIT15 = 0x8000;

const unsigned int CLR_BIT0 = 0xFFFFFFFE;
const unsigned int CLR_BIT1 = 0xFFFFFFFD;
const unsigned int CLR_BIT2 = 0xFFFFFFFB;
const unsigned int CLR_BIT3 = 0xFFFFFFF7;
const unsigned int CLR_BIT4 = 0xFFFFFFEF;
const unsigned int CLR_BIT5 = 0xFFFFFFDF;
const unsigned int CLR_BIT6 = 0xFFFFFFBF;
const unsigned int CLR_BIT7 = 0xFFFFFF7F;
const unsigned int CLR_BIT8 = 0xFFFFFEFF;
const unsigned int CLR_BIT9 = 0xFFFFFDFF;
const unsigned int CLR_BIT10 = 0xFFFFFBFF;
const unsigned int CLR_BIT11 = 0xFFFFF7FF;
const unsigned int CLR_BIT12 = 0xFFFFEFFF;
const unsigned int CLR_BIT13 = 0xFFFFDFFF;
const unsigned int CLR_BIT14 = 0xFFFFBFFF;
const unsigned int CLR_BIT15 = 0xFFFF7FFF;
const unsigned int CLR_BIT16 = 0xFFFEFFFF;
const unsigned int CLR_BIT17 = 0xFFFDFFFF;
const unsigned int CLR_BIT18 = 0xFFFBFFFF;
const unsigned int CLR_BIT19 = 0xFFF7FFFF;
const unsigned int CLR_BIT20 = 0xFFEFFFFF;
const unsigned int CLR_BIT21 = 0xFFDFFFFF;
const unsigned int CLR_BIT22 = 0xFFBFFFFF;
const unsigned int CLR_BIT23 = 0xFF7FFFFF;
const unsigned int CLR_BIT24 = 0xFEFFFFFF;
const unsigned int CLR_BIT25 = 0xFDFFFFFF;
const unsigned int CLR_BIT26 = 0xFBFFFFFF;
const unsigned int CLR_BIT27 = 0xF7FFFFFF;
const unsigned int CLR_BIT28 = 0xEFFFFFFF;
const unsigned int CLR_BIT29 = 0xDFFFFFFF;
const unsigned int CLR_BIT30 = 0xBFFFFFFF;
const unsigned int CLR_BIT31 = 0x7FFFFFFF;

const unsigned int SET_BIT0 = 0x00000001;
const unsigned int SET_BIT1 = 0x00000002;
const unsigned int SET_BIT2 = 0x00000004;
const unsigned int SET_BIT3 = 0x00000008;
const unsigned int SET_BIT4 = 0x00000010;
const unsigned int SET_BIT5 = 0x00000020;
const unsigned int SET_BIT6 = 0x00000040;
const unsigned int SET_BIT7 = 0x00000080;
const unsigned int SET_BIT8 = 0x00000100;
const unsigned int SET_BIT9 = 0x00000200;
const unsigned int SET_BIT10 = 0x00000400;
const unsigned int SET_BIT11 = 0x00000800;
const unsigned int SET_BIT12 = 0x00001000;
const unsigned int SET_BIT13 = 0x00002000;
const unsigned int SET_BIT14 = 0x00004000;
const unsigned int SET_BIT15 = 0x00008000;
const unsigned int SET_BIT16 = 0x00010000;
const unsigned int SET_BIT17 = 0x00020000;
const unsigned int SET_BIT18 = 0x00040000;
const unsigned int SET_BIT19 = 0x00080000;
const unsigned int SET_BIT20 = 0x00100000;
const unsigned int SET_BIT21 = 0x00200000;
const unsigned int SET_BIT22 = 0x00400000;
const unsigned int SET_BIT23 = 0x00800000;
const unsigned int SET_BIT24 = 0x01000000;
const unsigned int SET_BIT25 = 0x02000000;
const unsigned int SET_BIT26 = 0x04000000;
const unsigned int SET_BIT27 = 0x08000000;
const unsigned int SET_BIT28 = 0x10000000;
const unsigned int SET_BIT29 = 0x20000000;
const unsigned int SET_BIT30 = 0x40000000;
const unsigned int SET_BIT31 = 0x80000000;

const unsigned int BIT0_MASK = 0x00000001;
const unsigned int BIT1_MASK = 0x00000002;
const unsigned int BIT2_MASK = 0x00000004;
const unsigned int BIT3_MASK = 0x00000008;
const unsigned int BIT4_MASK = 0x00000010;
const unsigned int BIT5_MASK = 0x00000020;
const unsigned int BIT6_MASK = 0x00000040;
const unsigned int BIT7_MASK = 0x00000080;
const unsigned int BIT8_MASK = 0x00000100;
const unsigned int BIT9_MASK = 0x00000200;
const unsigned int BIT10_MASK = 0x00000400;
const unsigned int BIT11_MASK = 0x00000800;
const unsigned int BIT12_MASK = 0x00001000;
const unsigned int BIT13_MASK = 0x00002000;
const unsigned int BIT14_MASK = 0x00004000;
const unsigned int BIT15_MASK = 0x00008000;
const unsigned int BIT16_MASK = 0x00010000;
const unsigned int BIT17_MASK = 0x00020000;
const unsigned int BIT18_MASK = 0x00040000;
const unsigned int BIT19_MASK = 0x00080000;
const unsigned int BIT20_MASK = 0x00100000;
const unsigned int BIT21_MASK = 0x00200000;
const unsigned int BIT22_MASK = 0x00400000;
const unsigned int BIT23_MASK = 0x00800000;
const unsigned int BIT24_MASK = 0x01000000;
const unsigned int BIT25_MASK = 0x02000000;
const unsigned int BIT26_MASK = 0x04000000;
const unsigned int BIT27_MASK = 0x08000000;
const unsigned int BIT28_MASK = 0x10000000;
const unsigned int BIT29_MASK = 0x20000000;
const unsigned int BIT30_MASK = 0x40000000;
const unsigned int BIT31_MASK = 0x80000000;

//===========================================================
//  pseudo-functions (JTA)
//===========================================================

#define EXTRACT_BIT0(x) (x & BIT0_MASK)
#define EXTRACT_BIT1(x) ((x & BIT1_MASK) >> 1)
#define EXTRACT_BIT2(x) ((x & BIT2_MASK) >> 2)
#define EXTRACT_BIT3(x) ((x & BIT3_MASK) >> 3)
#define EXTRACT_BIT4(x) ((x & BIT4_MASK) >> 4)
#define EXTRACT_BIT5(x) ((x & BIT5_MASK) >> 5)
#define EXTRACT_BIT6(x) ((x & BIT6_MASK) >> 6)
#define EXTRACT_BIT7(x) ((x & BIT7_MASK) >> 7)
#define EXTRACT_BIT8(x) ((x & BIT8_MASK) >> 8)
#define EXTRACT_BIT9(x) ((x & BIT9_MASK) >> 9)
#define EXTRACT_BIT10(x) ((x & BIT10_MASK) >> 10)
#define EXTRACT_BIT11(x) ((x & BIT11_MASK) >> 11)
#define EXTRACT_BIT12(x) ((x & BIT12_MASK) >> 12)
#define EXTRACT_BIT13(x) ((x & BIT13_MASK) >> 13)
#define EXTRACT_BIT14(x) ((x & BIT14_MASK) >> 14)
#define EXTRACT_BIT15(x) ((x & BIT15_MASK) >> 15)
#define EXTRACT_BIT16(x) ((x & BIT16_MASK) >> 16)
#define EXTRACT_BIT17(x) ((x & BIT17_MASK) >> 17)
#define EXTRACT_BIT18(x) ((x & BIT18_MASK) >> 18)
#define EXTRACT_BIT19(x) ((x & BIT19_MASK) >> 19)
#define EXTRACT_BIT20(x) ((x & BIT20_MASK) >> 20)
#define EXTRACT_BIT21(x) ((x & BIT21_MASK) >> 21)
#define EXTRACT_BIT22(x) ((x & BIT22_MASK) >> 22)
#define EXTRACT_BIT23(x) ((x & BIT23_MASK) >> 23)
#define EXTRACT_BIT24(x) ((x & BIT24_MASK) >> 24)
#define EXTRACT_BIT25(x) ((x & BIT25_MASK) >> 25)
#define EXTRACT_BIT26(x) ((x & BIT26_MASK) >> 26)
#define EXTRACT_BIT27(x) ((x & BIT27_MASK) >> 27)
#define EXTRACT_BIT28(x) ((x & BIT28_MASK) >> 28)
#define EXTRACT_BIT29(x) ((x & BIT29_MASK) >> 22)
#define EXTRACT_BIT30(x) ((x & BIT30_MASK) >> 30)
#define EXTRACT_BIT31(x) ((x & BIT31_MASK) >> 31)

#define READ_MOD_WRITE_BIT0(x)  CLR_BIT0, ((x & 0x1) << 0)
#define READ_MOD_WRITE_BIT1(x)  CLR_BIT1, ((x & 0x1) << 1)
#define READ_MOD_WRITE_BIT2(x)  CLR_BIT2, ((x & 0x1) << 2)
#define READ_MOD_WRITE_BIT3(x)  CLR_BIT3, ((x & 0x1) << 3)
#define READ_MOD_WRITE_BIT4(x)  CLR_BIT4, ((x & 0x1) << 4)
#define READ_MOD_WRITE_BIT5(x)  CLR_BIT5, ((x & 0x1) << 5)
#define READ_MOD_WRITE_BIT6(x)  CLR_BIT6, ((x & 0x1) << 6)
#define READ_MOD_WRITE_BIT7(x)  CLR_BIT7, ((x & 0x1) << 7)
#define READ_MOD_WRITE_BIT8(x)  CLR_BIT8, ((x & 0x1) << 8)
#define READ_MOD_WRITE_BIT9(x)  CLR_BIT9, ((x & 0x1) << 9)
#define READ_MOD_WRITE_BIT10(x) CLR_BIT10,((x & 0x1) << 10)
#define READ_MOD_WRITE_BIT11(x) CLR_BIT11,((x & 0x1) << 11)
#define READ_MOD_WRITE_BIT12(x) CLR_BIT12,((x & 0x1) << 12)
#define READ_MOD_WRITE_BIT13(x) CLR_BIT13,((x & 0x1) << 13)
#define READ_MOD_WRITE_BIT14(x) CLR_BIT14,((x & 0x1) << 14)
#define READ_MOD_WRITE_BIT15(x) CLR_BIT15,((x & 0x1) << 15)
#define READ_MOD_WRITE_BIT16(x) CLR_BIT16,((x & 0x1) << 16)
#define READ_MOD_WRITE_BIT17(x) CLR_BIT17,((x & 0x1) << 17)
#define READ_MOD_WRITE_BIT18(x) CLR_BIT18,((x & 0x1) << 18)
#define READ_MOD_WRITE_BIT19(x) CLR_BIT19,((x & 0x1) << 19)
#define READ_MOD_WRITE_BIT20(x) CLR_BIT20,((x & 0x1) << 20)
#define READ_MOD_WRITE_BIT21(x) CLR_BIT21,((x & 0x1) << 21)
#define READ_MOD_WRITE_BIT22(x) CLR_BIT22,((x & 0x1) << 22)
#define READ_MOD_WRITE_BIT23(x) CLR_BIT23,((x & 0x1) << 23)
#define READ_MOD_WRITE_BIT24(x) CLR_BIT24,((x & 0x1) << 24)
#define READ_MOD_WRITE_BIT25(x) CLR_BIT25,((x & 0x1) << 25)
#define READ_MOD_WRITE_BIT26(x) CLR_BIT26,((x & 0x1) << 26)
#define READ_MOD_WRITE_BIT27(x) CLR_BIT27,((x & 0x1) << 27)
#define READ_MOD_WRITE_BIT28(x) CLR_BIT28,((x & 0x1) << 28)
#define READ_MOD_WRITE_BIT29(x) CLR_BIT29,((x & 0x1) << 29)
#define READ_MOD_WRITE_BIT30(x) CLR_BIT30,((x & 0x1) << 30)
#define READ_MOD_WRITE_BIT31(x) CLR_BIT31,((x & 0x1) << 31)

extern PARS Pars;
extern long long int dgsHeaderID[20];


/*------------------------------------------------------------------*/

int
DGSEvDecompose_v3 (unsigned int *ev, int len, DGSEVENT * DGSEvent, int tlkup[], int tid[])
{

  /* declarations */

  int i, k;
  unsigned int ui0 = 0;
  unsigned int t1 = 0, t2 = 0, t3 = 0, t4 = 0;
  unsigned long long int ulli1;


  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("entered DGSEvDecompose_v3:\n");

  /* swap the bytes */

  i = 0;
  while (i < len)
    {

      /* before 4 3 2 1 */
      /*        | | | | */
      /* after  1 2 3 4 */

      t1 = (*(ev + i) & 0x000000ff) << 24;
      t2 = (*(ev + i) & 0x0000ff00) << 8;
      t3 = (*(ev + i) & 0x00ff0000) >> 8;
      t4 = (*(ev + i) & 0xff000000) >> 24;
      *(ev + i) = t1 + t2 + t3 + t4;

      i++;
    }

  /* debug print */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("event len=%i (%lui bytes) >\n", len, len * sizeof (unsigned int));
      for (i = 0; i < len; i++)
        {
          printf ("%3i[doc: %3i]: %12u, 0x%8.8x; ", i, i + 1, *(ev + i), *(ev + i));
          ui0 = 0x80000000;
          printf ("|");
          for (k = 0; k <= 31; k++)
            {
              if ((*(ev + i) & ui0) == ui0)
                printf ("1");
              else
                printf ("0");
              ui0 = ui0 / 2;
              if ((k + 1) % 4 == 0)
                printf ("|");
            };
          printf ("\n");
        };
    };

  // Decode the generic part of the header.

  DGSEvent->chan_id = (*(ev + 0) & 0x0000000f);
  DGSEvent->board_id = ((*(ev + 0) & 0x0000fff0) >> 4); // USER_DEF
  DGSEvent->id = DGSEvent->board_id * 10 + DGSEvent->chan_id;
  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("chan_id= %i - ", DGSEvent->chan_id);
      printf ("board_id= %i - ", DGSEvent->board_id);
      printf ("id= %i\n", DGSEvent->id);
      printf ("\n");
    };
  DGSEvent->packet_length = ((*(ev + 0) & 0x07ff0000) >> 16);
  DGSEvent->geo_addr = ((*(ev + 0) & 0xf8000000) >> 27);

  DGSEvent->header_type = ((*(ev + 2) & 0x000f0000) >> 16);
  dgsHeaderID[DGSEvent->header_type]++;
  DGSEvent->event_type = ((*(ev + 2) & 0x03800000) >> 23);      // hope this is right.
  DGSEvent->header_length = ((*(ev + 2) & 0xFC000000) >> 26);   // hope this is right.

  /* extract the LED time stamp */

  DGSEvent->event_timestamp = 0;
  DGSEvent->event_timestamp = (unsigned long long int) *(ev + 1);
  ulli1 = (unsigned long long int) (*(ev + 2) & 0x0000ffff);
  ulli1 = (ulli1 << 32);
  DGSEvent->event_timestamp += ulli1;

  /* store the type and type ID */

  DGSEvent->tpe = tlkup[DGSEvent->id];
  DGSEvent->tid = tid[DGSEvent->id];
  DGSEvent->flag = 0;
//  for (i=1000;i<1030;i++)
//    printf("debug2: %i: tlkup %i tid %i\n", i, tlkup[i], tid[i]);

  /* extract the header values based on header type */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("DGSEvent->header_type=%i\n", DGSEvent->header_type);
  switch (DGSEvent->header_type)
    {

    case 0:                    // Original LED header (Compatibility mode)
      DGSEvent->timestamp_match_flag = 0;
      DGSEvent->external_disc_flag = 0;
      DGSEvent->cfd_valid_flag = 0;
      DGSEvent->pileup_only_flag = 0;
      DGSEvent->cfd_sample_0 = 0;
      DGSEvent->cfd_sample_1 = 0;
      DGSEvent->cfd_sample_2 = 0;

      DGSEvent->peak_valid_flag = ((*(ev + 3) & 0x00000200) >> 9);      // Word 3: 9
      DGSEvent->offset_flag = ((*(ev + 3) & 0x00000400) >> 10); // Word 3: 10
      DGSEvent->sync_error_flag = ((*(ev + 3) & 0x00001000) >> 12);     // Word 3: 12
      DGSEvent->general_error_flag = ((*(ev + 3) & 0x00002000) >> 13);  // Word 3: 13
      DGSEvent->pileup_flag = ((*(ev + 3) & 0x00008000) >> 15); // Word 3: 15
      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 & 
        (((unsigned long long int) (*(ev + 4) & 0xFFFFFFFF)) << 16);    // Word 4 :31..0  
      DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 0);     // Word 5: 23..0
      DGSEvent->sum1 = ((*(ev + 7) & 0x00FFFFFF) >> 0); // Word 7: 23..0
      DGSEvent->sum2 = ((*(ev + 7) & 0xFF000000) >> 28) |       // Word 7: 31..24 & 
        ((*(ev + 8) & 0x0000FFFF) << 8);        // Word 8: 15..0 
      DGSEvent->peak_timestamp = (((unsigned long long int) (*(ev + 8) & 0xFFFF0000)) >> 16) |  // Word 8: 31..16 & 
        (((unsigned long long int) (*(ev + 9) & 0xFFFFFFFF)) << 16);    // Word 9 :31..0  
      DGSEvent->m2begin = ((*(ev + 10) & 0x00003FFF) >> 0);     // Word 10:13..0
      DGSEvent->m2end = ((*(ev + 10) & 0x3FFF0000) >> 16);      // Word 10:29..16
      DGSEvent->m1begin = ((*(ev + 11) & 0x00003FFF) >> 0);     // Word 11:13..0
      DGSEvent->m1end = ((*(ev + 11) & 0x3FFF0000) >> 16);      // Word 11:29..16
      DGSEvent->peak_sample = ((*(ev + 12) & 0x00003FFF) >> 0); // Word 12:13..0
      DGSEvent->base_sample = ((*(ev + 12) & 0x3FFF0000) >> 16);        // Word 12:29..16
      break;

    case 1:                    // New LED Header
      DGSEvent->timestamp_match_flag = 0;
      DGSEvent->cfd_valid_flag = 0;
      DGSEvent->cfd_sample_0 = 0;
      DGSEvent->cfd_sample_1 = 0;
      DGSEvent->cfd_sample_2 = 0;

      DGSEvent->external_disc_flag = ((*(ev + 3) & 0x00000100) >> 8);   // Word 3: 8
      DGSEvent->peak_valid_flag = ((*(ev + 3) & 0x00000200) >> 9);      // Word 3: 9
      DGSEvent->offset_flag = ((*(ev + 3) & 0x00000400) >> 10); // Word 3: 10
      DGSEvent->sync_error_flag = ((*(ev + 3) & 0x00001000) >> 12);     // Word 3: 12
      DGSEvent->general_error_flag = ((*(ev + 3) & 0x00002000) >> 13);  // Word 3: 13
      DGSEvent->pileup_only_flag = ((*(ev + 3) & 0x00004000) >> 14);    // Word 3: 14
      DGSEvent->pileup_flag = ((*(ev + 3) & 0x00008000) >> 15); // Word 3: 15  
      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 & 
        (((unsigned long long int) (*(ev + 4) & 0xFFFFFFFF)) << 16);    // Word 4 :31..0  
      DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 0);     // Word 5: 23..0
      DGSEvent->sum1 = ((*(ev + 7) & 0x00FFFFFF) >> 0); // Word 7: 23..0
      DGSEvent->sum2 = ((*(ev + 7) & 0xFF000000) >> 28) |       // Word 7: 31..24 & 
        ((*(ev + 8) & 0x0000FFFF) << 8);        // Word 8: 15..0 
      DGSEvent->peak_timestamp = (((unsigned long long int) (*(ev + 8) & 0xFFFF0000)) >> 16) |  // Word 8: 31..16 & 
        (((unsigned long long int) (*(ev + 9) & 0xFFFFFFFF)) << 16);    // Word 9: 31..0  
      DGSEvent->m1begin = ((*(ev + 10) & 0x00003FFF) >> 0);     // Word 10:13..0
      //DGSEvent->m1end          = ((*(ev + 10) & 0x3FFF0000) >> 16);                             // Word 10:29..16

      DGSEvent->m2last_end_sample = ((*(ev + 9) & 0x00003FFF) >> 0);    // Word 9:13..0  
      DGSEvent->m2last_begin_sample = ((*(ev + 10) & 0x3FFF0000) >> 16);        // Word 10:13..0

      DGSEvent->m2end = ((*(ev + 11) & 0x00003FFF) >> 0);       // Word 11:13..0
      DGSEvent->m2begin = ((*(ev + 11) & 0x3FFF0000) >> 16);    // Word 11:29..16
      DGSEvent->peak_sample = ((*(ev + 12) & 0x00003FFF) >> 0); // Word 12:13..0
      DGSEvent->base_sample = ((*(ev + 12) & 0x3FFF0000) >> 16);        // Word 12:29..16
      break;

    case 3:                    // New LED Header
      DGSEvent->timestamp_match_flag = 0;
      DGSEvent->cfd_valid_flag = 0;
      DGSEvent->cfd_sample_0 = 0;
      DGSEvent->cfd_sample_1 = 0;
      DGSEvent->cfd_sample_2 = 0;

      DGSEvent->external_disc_flag = ((*(ev + 3) & 0x00000100) >> 8);   // Word 3: 8
      DGSEvent->peak_valid_flag = ((*(ev + 3) & 0x00000200) >> 9);      // Word 3: 9
      DGSEvent->offset_flag = ((*(ev + 3) & 0x00000400) >> 10); // Word 3: 10
      DGSEvent->sync_error_flag = ((*(ev + 3) & 0x00001000) >> 12);     // Word 3: 12
      DGSEvent->general_error_flag = ((*(ev + 3) & 0x00002000) >> 13);  // Word 3: 13
      DGSEvent->pileup_only_flag = ((*(ev + 3) & 0x00004000) >> 14);    // Word 3: 14
      DGSEvent->pileup_flag = ((*(ev + 3) & 0x00008000) >> 15); // Word 3: 15  
      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 & 
        (((unsigned long long int) (*(ev + 4) & 0xFFFFFFFF)) << 16);    // Word 4 :31..0  
      DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 0);     // Word 5: 23..0
      DGSEvent->sum1 = ((*(ev + 7) & 0x00FFFFFF) >> 0); // Word 7: 23..0
      DGSEvent->sum2 = ((*(ev + 7) & 0xFF000000) >> 28) |       // Word 7: 31..24 & 
        ((*(ev + 8) & 0x0000FFFF) << 8);        // Word 8: 15..0 
      DGSEvent->peak_timestamp = (((unsigned long long int) (*(ev + 8) & 0xFFFF0000)) >> 16) |  // Word 8: 31..16 & 
        (((unsigned long long int) (*(ev + 9) & 0xFFFFFFFF)) << 16);    // Word 9: 31..0  
      DGSEvent->m1begin = ((*(ev + 10) & 0x00003FFF) >> 0);     // Word 10:13..0
      //DGSEvent->m1end          = ((*(ev + 10) & 0x3FFF0000) >> 16);                             // Word 10:29..16

      DGSEvent->m2last_end_sample = ((*(ev + 9) & 0x00003FFF) >> 0);    // Word 9:13..0  
      DGSEvent->m2last_begin_sample = ((*(ev + 10) & 0x3FFF0000) >> 16);        // Word 10:13..0

      DGSEvent->m2end = ((*(ev + 11) & 0x00003FFF) >> 0);       // Word 11:13..0
      DGSEvent->m2begin = ((*(ev + 11) & 0x3FFF0000) >> 16);    // Word 11:29..16
      DGSEvent->peak_sample = ((*(ev + 12) & 0x00003FFF) >> 0); // Word 12:13..0
      DGSEvent->base_sample = ((*(ev + 12) & 0x3FFF0000) >> 16);        // Word 12:29..16
      break;

    case 4:                    // was 2              // CFD Header
      DGSEvent->timestamp_match_flag = ((*(ev + 3) & 0x00000080) >> 7); // Word 3: 7
      DGSEvent->external_disc_flag = ((*(ev + 3) & 0x00000100) >> 8);   // Word 3: 8
      DGSEvent->peak_valid_flag = ((*(ev + 3) & 0x00000200) >> 9);      // Word 3: 9
      DGSEvent->offset_flag = ((*(ev + 3) & 0x00000400) >> 10); // Word 3: 10
      DGSEvent->cfd_valid_flag = ((*(ev + 3) & 0x00000800) >> 11);      // Word 3: 11
      DGSEvent->sync_error_flag = ((*(ev + 3) & 0x00001000) >> 12);     // Word 3: 12
      DGSEvent->general_error_flag = ((*(ev + 3) & 0x00002000) >> 13);  // Word 3: 13
      DGSEvent->pileup_only_flag = ((*(ev + 3) & 0x00004000) >> 14);    // Word 3: 14
      DGSEvent->pileup_flag = ((*(ev + 3) & 0x00008000) >> 15); // Word 3: 15
      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 &
        (((unsigned long long int) (*(ev + 4) & 0xFFFFFFFF)) << 16);    // Word 4 :31..0
      DGSEvent->cfd_sample_0 = ((*(ev + 4) & 0x3FFF0000) >> 16);        // Word 4: 29..16
      DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 0);     // Word 5: 23..0
      DGSEvent->cfd_sample_1 = ((*(ev + 6) & 0x00003FFF) >> 0); // Word 6: 13..0
      DGSEvent->cfd_sample_2 = ((*(ev + 6) & 0x3FFF0000) >> 16);        // Word 6: 29..16
      DGSEvent->sum1 = ((*(ev + 7) & 0x00FFFFFF) >> 0); // Word 7: 23..0
      DGSEvent->sum2 = ((*(ev + 7) & 0xFF000000) >> 24) |       //was 28                      // Word 7: 31..24 & 
        ((*(ev + 8) & 0x0000FFFF) << 8);        // Word 8: 15..0 
      DGSEvent->peak_timestamp = (((unsigned long long int) (*(ev + 8) & 0xFFFF0000)) >> 16) |  // Word 8: 31..16 & 
        (((unsigned long long int) (*(ev + 9) & 0xFFFFFFFF)) << 16);    // Word 9: 31..0  

      DGSEvent->m2end = ((*(ev + 10) & 0x00003FFF) >> 0);       // Word 10:13..0   POST_RISE_SAMPLE
      DGSEvent->m2begin = ((*(ev + 10) & 0x3FFF0000) >> 16);    // Word 10:29..16  LAST_POST_RISE_SAMPLE
      DGSEvent->m1begin = ((*(ev + 11) & 0x00003FFF) >> 0);     // Word 11:13..0   PRE_RISE_ENTER_SAMPLE
      DGSEvent->m1end = ((*(ev + 11) & 0x3FFF0000) >> 16);      // Word 11:29..16  PRE_RISE_LEAVE_SAMPLE

      DGSEvent->peak_sample = ((*(ev + 12) & 0x00003FFF) >> 0); // Word 12:13..0
      DGSEvent->base_sample = ((*(ev + 12) & 0x3FFF0000) >> 16);        // Word 12:29..16
      break;

    case 5:                    // New LED Header //TL: copy of 1 for now
      DGSEvent->timestamp_match_flag = 0;
      DGSEvent->cfd_valid_flag = 0;
      DGSEvent->cfd_sample_0 = 0;
      DGSEvent->cfd_sample_1 = 0;
      DGSEvent->cfd_sample_2 = 0;

      DGSEvent->external_disc_flag = ((*(ev + 3) & 0x00000100) >> 8);   // Word 3: 8
      DGSEvent->peak_valid_flag = ((*(ev + 3) & 0x00000200) >> 9);      // Word 3: 9
      DGSEvent->offset_flag = ((*(ev + 3) & 0x00000400) >> 10); // Word 3: 10
      DGSEvent->sync_error_flag = ((*(ev + 3) & 0x00001000) >> 12);     // Word 3: 12
      DGSEvent->general_error_flag = ((*(ev + 3) & 0x00002000) >> 13);  // Word 3: 13
      DGSEvent->pileup_only_flag = ((*(ev + 3) & 0x00004000) >> 14);    // Word 3: 14
      DGSEvent->pileup_flag = ((*(ev + 3) & 0x00008000) >> 15); // Word 3: 15  
      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 & 
        (((unsigned long long int) (*(ev + 4) & 0xFFFFFFFF)) << 16);    // Word 4 :31..0  
      DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 0);     // Word 5: 23..0
      DGSEvent->sum1 = ((*(ev + 7) & 0x00FFFFFF) >> 0); // Word 7: 23..0
      DGSEvent->sum2 = ((*(ev + 7) & 0xFF000000) >> 28) |       // Word 7: 31..24 & 
        ((*(ev + 8) & 0x0000FFFF) << 8);        // Word 8: 15..0 
      DGSEvent->peak_timestamp = (((unsigned long long int) (*(ev + 8) & 0xFFFF0000)) >> 16) |  // Word 8: 31..16 & 
        (((unsigned long long int) (*(ev + 9) & 0xFFFFFFFF)) << 16);    // Word 9: 31..0  
      DGSEvent->m1begin = ((*(ev + 10) & 0x00003FFF) >> 0);     // Word 10:13..0
      //DGSEvent->m1end          = ((*(ev + 10) & 0x3FFF0000) >> 16);                             // Word 10:29..16

      DGSEvent->m2last_end_sample = ((*(ev + 9) & 0x00003FFF) >> 0);    // Word 9:13..0  
      DGSEvent->m2last_begin_sample = ((*(ev + 10) & 0x3FFF0000) >> 16);        // Word 10:13..0

      DGSEvent->m2end = ((*(ev + 11) & 0x00003FFF) >> 0);       // Word 11:13..0
      DGSEvent->m2begin = ((*(ev + 11) & 0x3FFF0000) >> 16);    // Word 11:29..16
      DGSEvent->peak_sample = ((*(ev + 12) & 0x00003FFF) >> 0); // Word 12:13..0
      DGSEvent->base_sample = ((*(ev + 12) & 0x3FFF0000) >> 16);        // Word 12:29..16
      break;

    case 6:                    // was 2              // CFD Header   //TL: copy of 2 for now
      DGSEvent->timestamp_match_flag = ((*(ev + 3) & 0x00000080) >> 7); // Word 3: 7
      DGSEvent->external_disc_flag = ((*(ev + 3) & 0x00000100) >> 8);   // Word 3: 8
      DGSEvent->peak_valid_flag = ((*(ev + 3) & 0x00000200) >> 9);      // Word 3: 9
      DGSEvent->offset_flag = ((*(ev + 3) & 0x00000400) >> 10); // Word 3: 10
      DGSEvent->cfd_valid_flag = ((*(ev + 3) & 0x00000800) >> 11);      // Word 3: 11
      DGSEvent->sync_error_flag = ((*(ev + 3) & 0x00001000) >> 12);     // Word 3: 12
      DGSEvent->general_error_flag = ((*(ev + 3) & 0x00002000) >> 13);  // Word 3: 13
      DGSEvent->pileup_only_flag = ((*(ev + 3) & 0x00004000) >> 14);    // Word 3: 14
      DGSEvent->pileup_flag = ((*(ev + 3) & 0x00008000) >> 15); // Word 3: 15
      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 &
        (((unsigned long long int) (*(ev + 4) & 0xFFFFFFFF)) << 16);    // Word 4 :31..0
      DGSEvent->cfd_sample_0 = ((*(ev + 4) & 0x3FFF0000) >> 16);        // Word 4: 29..16
      DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 0);     // Word 5: 23..0
      DGSEvent->cfd_sample_1 = ((*(ev + 6) & 0x00003FFF) >> 0); // Word 6: 13..0
      DGSEvent->cfd_sample_2 = ((*(ev + 6) & 0x3FFF0000) >> 16);        // Word 6: 29..16
      DGSEvent->sum1 = ((*(ev + 7) & 0x00FFFFFF) >> 0); // Word 7: 23..0
      DGSEvent->sum2 = ((*(ev + 7) & 0xFF000000) >> 24) |       //was 28                      // Word 7: 31..24 & 
        ((*(ev + 8) & 0x0000FFFF) << 8);        // Word 8: 15..0 
      DGSEvent->peak_timestamp = (((unsigned long long int) (*(ev + 8) & 0xFFFF0000)) >> 16) |  // Word 8: 31..16 & 
        (((unsigned long long int) (*(ev + 9) & 0xFFFFFFFF)) << 16);    // Word 9: 31..0  

      DGSEvent->m2end = ((*(ev + 10) & 0x00003FFF) >> 0);       // Word 10:13..0   POST_RISE_SAMPLE
      DGSEvent->m2begin = ((*(ev + 10) & 0x3FFF0000) >> 16);    // Word 10:29..16  LAST_POST_RISE_SAMPLE
      DGSEvent->m1begin = ((*(ev + 11) & 0x00003FFF) >> 0);     // Word 11:13..0   PRE_RISE_ENTER_SAMPLE
      DGSEvent->m1end = ((*(ev + 11) & 0x3FFF0000) >> 16);      // Word 11:29..16  PRE_RISE_LEAVE_SAMPLE

      DGSEvent->peak_sample = ((*(ev + 12) & 0x00003FFF) >> 0); // Word 12:13..0
      DGSEvent->base_sample = ((*(ev + 12) & 0x3FFF0000) >> 16);        // Word 12:29..16
      break;

    default:
      printf ("DGSEvent->header_type= %i is unknown to this version of bin_dgs.c; \n", DGSEvent->header_type);
      printf ("current ev no: %i; header dump follows:\n", Pars.CurEvNo);
      printf ("event len=%i (%lui bytes) >\n", len, len * sizeof (unsigned int));
      for (i = 0; i < len; i++)
        {
          printf ("%3i[doc: %3i]: %12u, 0x%8.8x; ", i, i + 1, *(ev + i), *(ev + i));
          ui0 = 0x80000000;
          printf ("|");
          for (k = 0; k <= 31; k++)
            {
              if ((*(ev + i) & ui0) == ui0)
                printf ("1");
              else
                printf ("0");
              ui0 = ui0 / 2;
              if ((k + 1) % 4 == 0)
                printf ("|");
            };
          printf ("\n");
        };
      return (1);
      break;

    }

  DGSEvent->baseline = (DGSEvent->m2end + DGSEvent->m2begin) / 2;


  /* Now load Trace into DGSEvent Structure */

  DGSEvent->traceLen = 0;
  for (i = 13; i < len - 1; i++)
    {
      if (i < 1037)
        {
          DGSEvent->trace[2 * (i - 13)] = (unsigned short int) (*(ev + i) & 0xffff);
          DGSEvent->trace[2 * (i - 13) + 1] = (unsigned short int) ((*(ev + i) >> 16) & 0xffff);
          DGSEvent->traceLen += 2;
        }
    }


  /* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit DGSEvDecompose_v3\n");

  return (0);

}

