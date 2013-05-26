/*------------------------------------------------------------------------------
* sirf.c : SiRF III/IV receiver dependent functions
*
*    Copyright (C) 2013 A.Illarionov, All rights reserved.
*
*    [1] SiRFstarIV One Socket Protocol Interface Control Document 
*
* version : $Revision:$ $Date:$
* history : 2013/05/26 1.0  first version
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SIRFSYNC1    0xA0        /* SiRF binary message start sync code 1 */
#define SIRFSYNC2    0xA2        /* SiRF binary message start sync code 2 */
#define SIRFSYNCEND1 0xB0        /* SiRF binary message end sync code 1 */
#define SIRFSYNCEND2 0xB3        /* SiRF binary message end sync code 2 */

#define MID_SRFCLOCK  0x07        /* SiRF binary clock status */
#define MID_SRF50BPS  0x08        /* SiRF binary 50BPS data */
#define MID_SRFNLMEAS 0x1c        /* SiRF binary navlib measurement data */

#define SIRFMAXRAWLEN 2047       /* max length of SiRF binary raw message */

/* get fields (big-endian) ---------------------------------------------------*/
#define U1(p) (*((unsigned char *)(p)))
#define I1(p) (*((char *)(p)))

static unsigned short U2(unsigned char *p)
{
    unsigned short value;
    unsigned char *q=(unsigned char *)&value+1;
    int i;
    for (i=0;i<2;i++) *q--=*p++;
    return value;
}
static unsigned int U4(unsigned char *p)
{
    unsigned int value;
    unsigned char *q=(unsigned char *)&value+3;
    int i;
    for (i=0;i<4;i++) *q--=*p++;
    return value;
}
static float R4(unsigned char *p)
{
    float value;
    unsigned char *q=(unsigned char *)&value+3;
    int i;
    for (i=0;i<4;i++) *q--=*p++;
    return value;
}
static double R8(unsigned char *p, unsigned gsw230)
{
    double value;
    unsigned char *q;
    int i;
    if (gsw230) {
    double value;
    unsigned char *q=(unsigned char *)&value+7;
    int i;
    for (i=0;i<8;i++) *q--=*p++;
    return value;
    }else {
    q = (unsigned char *)&value+3;
    for (i=0;i<4;i++) *q--=*p++;
    q = (unsigned char *)&value+7;
    for (i=0;i<4;i++) *q--=*p++;
    return value;
    }
}
/* checksum ------------------------------------------------------------------*/
static int chksum(const unsigned char *buff, int len)
{
    int i;
    unsigned short sum=0;

    if (len<8) return 0;
    for (i=4;i<len-4;i++) sum=0x7fff&(sum+buff[i]);
    return (sum>>8)==buff[len-4]&&(sum&0xFF)==buff[len-3];
}
/* decode id#20 navigation data (user) ---------------------------------------*/
static int decode_sirfclock(raw_t *raw)
{
    unsigned drift;
    double bias;
    int i,j,n=0;
    unsigned char *p;
    unsigned wn;
    double tow;
    obsd_t *src,*dst;
    gtime_t time0={0};

    trace(4,"decode_sirfclock: len=%d\n",raw->len);

    p=raw->buff+4;
    if (raw->len!=28) {
        trace(2,"SiRF mid#7 length error: len=%d\n",raw->len);
        return -1;
    }

    wn=U2(p+1);
    tow=U4(p+3)*0.01;
    raw->time=gpst2time(wn,tow);
    drift=U4(p+8);
    bias=U4(p+12)/1.0e9;

    for (i=0;i<raw->obuf.n&&i<MAXOBS;i++) {
        src=&raw->obuf.data[i];
        if (!satsys(src->sat,NULL)) continue;
        if (fabs(tow+bias-src->time.sec)>0.1) continue;
        dst=&raw->obs.data[n++];
        *dst=*src;
        dst->time=gpst2time(wn,src->time.sec-bias);
        if (dst->L[0]) dst->L[0]-=FREQ1*bias;
        dst->P[0]-=CLIGHT*bias;
        dst->D[0]-=-1.0*drift;
    }
    raw->obs.n=n;

    /* clear observation data buffer */
    for (i=0;i<MAXOBS;i++) {
        raw->obuf.data[i].time=time0;
        for (j=0;j<NFREQ+NEXOBS;j++) {
            raw->obuf.data[i].L[j]=raw->obuf.data[i].P[j]=0.0;
            raw->obuf.data[i].D[j]=0.0;
            raw->obuf.data[i].SNR[j]=raw->obuf.data[i].LLI[j]=0;
            raw->obuf.data[i].code[j]=CODE_NONE;
        }
    }

    return n>0?1:0;
}
/* decode id#23 measurement block --------------------------------------------*/
static int decode_sirfnlmeas(raw_t *raw)
{
    unsigned ch,flags,pherrcnt;
    int i;
    unsigned char *p;
    unsigned gsw230;
    obsd_t *obsd;

    trace(4,"decode_sirfnlmeas: len=%d\n",raw->len);

    if (raw->len!=64) {
        trace(2,"SiRF mid#28 length error: len=%d\n",raw->len);
        return -1;
    }

    p=raw->buff+4;
    ch=U1(p+1);
    if (ch>=MAXOBS) {
        trace(2,"SiRF mid#28 wrong channel: ch=%d\n",ch);
        return -1;
    }
    gsw230=strstr(raw->opt,"-GSW230")!=NULL;
    obsd=&raw->obuf.data[ch];
    if (ch>=raw->obuf.n) raw->obuf.n=ch+1;
    obsd->sat=satno(SYS_GPS,U1(p+6));
    obsd->time.sec=R8(p+7,gsw230); /* XXX */
    obsd->P[0]=R8(p+15, gsw230);
    obsd->D[0]=-1*R4(p+23)/(CLIGHT/FREQ1);
    obsd->L[0]=R8(p+27,gsw230)/(CLIGHT/FREQ1);
    obsd->SNR[0]=U1(p+38);
    for (i=39;i<=47;i++) {
        if (U1(p+i)<obsd->SNR[0]) obsd->SNR[0]=U1(p+i);
    }
    obsd->SNR[0]*=4.0;
    flags=U1(p+37);
    pherrcnt=U1(p+54);
    obsd->LLI[0]=flags&0x02&&pherrcnt<50?0:1;
    obsd->code[0]=CODE_L1C;

    return 0;
}
/* save subframe -------------------------------------------------------------*/
static int save_subfrm(int sat, raw_t *raw)
{
    unsigned int i,id,word;
    unsigned char *p=raw->buff+7;
    unsigned char subfrm[30];

    for (i=0;i<10;i++) {
        word=U4(p+i*4);
        if (!decode_word(word,subfrm+i*3)) {
            trace(4,"save_subfrm: parity error sat=%2d i=%d word=%08x\n",sat,i,word);
            return 0;
        }
    }

    id=getbitu(subfrm,43,3);
    if (id<1||5<id) return 0;

    trace(4,"save_subfrm: sat=%2d id=%d\n",sat,id);
    memcpy(raw->subfrm[sat-1]+(id-1)*30,subfrm,30);

    return id;
}
/* decode ephemeris ----------------------------------------------------------*/
static int decode_ephem(int sat, raw_t *raw)
{
    eph_t eph={0};

    trace(4,"decode_ephem: sat=%2d\n",sat);

    if (decode_frame(raw->subfrm[sat-1]   ,&eph,NULL,NULL,NULL,NULL)!=1||
            decode_frame(raw->subfrm[sat-1]+30,&eph,NULL,NULL,NULL,NULL)!=2||
            decode_frame(raw->subfrm[sat-1]+60,&eph,NULL,NULL,NULL,NULL)!=3) return 0;

    if (!strstr(raw->opt,"-EPHALL")) {
        if (eph.iode==raw->nav.eph[sat-1].iode) return 0; /* unchanged */
    }
    eph.sat=sat;
    raw->nav.eph[sat-1]=eph;
    raw->ephsat=sat;
    return 2;
}
/* decode almanac and ion/utc ------------------------------------------------*/
static int decode_alm1(int sat, raw_t *raw)
{
    trace(4,"decode_alm1 : sat=%2d\n",sat);
    decode_frame(raw->subfrm[sat-1]+90,NULL,raw->nav.alm,raw->nav.ion_gps,
            raw->nav.utc_gps,&raw->nav.leaps);
    return 0;
}
/* decode almanac ------------------------------------------------------------*/
static int decode_alm2(int sat, raw_t *raw)
{
    trace(4,"decode_alm2 : sat=%2d\n",sat);
    decode_frame(raw->subfrm[sat-1]+120,NULL,raw->nav.alm,NULL,NULL,NULL);
    return  0;
}
/* decode SiRF 50BPS message ------------ ------------------------------------*/
static int decode_sirf50bps(raw_t *raw)
{
    int prn,sat,id;
    unsigned char *p=raw->buff+4;

    trace(4,"decode_sirf50bps: len=%d\n",raw->len);

    if (raw->len<51) {
        trace(2,"SiRF 50BSP length error: len=%d\n",raw->len);
        return -1;
    }
    prn=U1(p+2);
    if (!(sat=satno(SYS_GPS,prn))) {
        trace(2,"SiRF 50BPS satellite number error: prn=%d\n",prn);
        return -1;
    }

    id=save_subfrm(sat,raw);
    if (id==3) return decode_ephem(sat,raw);
    if (id==4) return decode_alm1 (sat,raw);
    if (id==5) return decode_alm2 (sat,raw);
    return 0;
}
/* decode SiRF raw message --------------------------------------------*/
static int decode_sirf(raw_t *raw)
{
    unsigned char *p=raw->buff;
    int mid=U1(p+4);

    trace(3,"decode_sirf: mid=%2d\n",mid);

    if (raw->buff[raw->len-2]!=SIRFSYNCEND1
            ||raw->buff[raw->len-1]!=SIRFSYNCEND2) {
        trace(2,"SiRF sync error");
        return -1;
    }
    if (!chksum(raw->buff,raw->len)) {
        trace(2,"SiRF message checksum error: mid=%d len=%d\n",mid,raw->len);
        return -1;
    }
    if (raw->outtype) {
        sprintf(raw->msgtype,"SiRF %2d (%4d):",mid,raw->len);
    }
    switch (mid) {
        case MID_SRFCLOCK:  return decode_sirfclock(raw);
        case MID_SRF50BPS:  return decode_sirf50bps(raw);
        case MID_SRFNLMEAS: return decode_sirfnlmeas(raw);
    }
    return 0;
}
/* sync code -----------------------------------------------------------------*/
static int sync_sirf(unsigned char *buff, unsigned char data)
{
    buff[0]=buff[1]; buff[1]=data;
    return buff[0]==SIRFSYNC1&&buff[1]==SIRFSYNC2;
}
/* input SiRF Star raw message from stream -----------------------------------
 * input next SiRF Star raw message from stream
 * args   : raw_t *raw   IO     receiver raw data control struct
 *          unsigned char data I stream data (1 byte)
 * return : status (-1: error message, 0: no message, 1: input observation data,
 *                  2: input ephemeris,
 *                  9: input ion/utc parameter)
 *-----------------------------------------------------------------------------*/
extern int input_sirf(raw_t *raw, unsigned char data)
{
    trace(5,"input_sirf: data=%02x\n",data);

    /* synchronize frame */
    if (raw->nbyte==0) {
        if (!sync_sirf(raw->buff,data)) return 0;
        raw->nbyte=2;
        return 0;
    }
    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte<4) return 0;

    if (raw->nbyte==4) {
        if ((raw->len=U2(raw->buff+2)+8)>SIRFMAXRAWLEN) {
            trace(2,"SiRF length error: len=%d\n",raw->len);
            raw->nbyte=0;
            return -1;
        }
    }

    if (raw->nbyte<4||raw->nbyte<raw->len) return 0;

    raw->nbyte=0;

    /* decode sirf raw message */
    return decode_sirf(raw);
}
/* input superstar 2 raw message from file -------------------------------------
 * input next superstar 2 raw message from file
 * args   : raw_t  *raw   IO     receiver raw data control struct
 *          FILE   *fp    I      file pointer
 * return : status(-2: end of file, -1...9: same as above)
 *-----------------------------------------------------------------------------*/
extern int input_sirff(raw_t *raw, FILE *fp)
{
    int i,data;

    trace(4,"input_sirff:\n");

    /* synchronize frame */
    if (raw->nbyte==0) {
        for (i=0;;i++) {
            if ((data=fgetc(fp))==EOF) return -2;
            if (sync_sirf(raw->buff,(unsigned char)data)) break;
            if (i>=4096) return 0;
        }
    }
    if (fread(raw->buff+2,1,2,fp)<2) return -2;
    raw->nbyte=4;

    if ((raw->len=U2(raw->buff+2)+8)>SIRFMAXRAWLEN) {
        trace(2,"SiRF length error: len=%d\n",raw->len);
        raw->nbyte=0;
        return -1;
    }
    if (fread(raw->buff+4,1,raw->len-4,fp)<(size_t)(raw->len-4)) return -2;
    raw->nbyte=0;

    /* decode SiRF raw message */
    return decode_sirf(raw);
}

