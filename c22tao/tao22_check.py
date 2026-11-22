#!/usr/bin/env python3
"""Replay round taoreplays22: closed-form bounds on pages of teorth/optimizationproblems (main 87b4ebb).
Each stated bound is evaluated by TWO separately written programs and compared with the printed digits.
 A: decimal arithmetic, 70 digits; pi by Machin, ln by atanh series.
 B: integer fixed point (10^70), pi by Gauss-Legendre, ln by AGM, sqrt by isqrt.
PASS needs |value - printed| <= 1 unit of the last printed digit (covers rounding and truncation) and A = B to 40 digits.
Red control: the printed value moved by 2 units must FAIL. Source: constant pages at teorth/optimizationproblems."""
from decimal import Decimal as D, getcontext
from math import isqrt
import sys
getcontext().prec=300
# ---------- program A ----------
def A_pi():
    def atan_inv(n):
        x=D(1)/n; s=x; t=x; k=1; n2=n*n
        while abs(t)>D(10)**-85:
            t=-t/n2; k+=2; s+=t/k
        return s
    return 4*(4*atan_inv(5)-atan_inv(239))
def A_ln(x):
    x=D(x); m=0
    while x>2: x/=2; m+=1
    while x<1: x*=2; m-=1
    def atanh_ln(y):  # ln y = 2 atanh((y-1)/(y+1))
        z=(y-1)/(y+1); z2=z*z; s=z; t=z; k=1
        while abs(t)>D(10)**-85:
            t*=z2; k+=2; s+=t/k
        return 2*s
    return atanh_ln(x/ D(1) ) + m*atanh_ln(D(2)) if m else atanh_ln(x)
PI=A_pi()
def A_val(i):
    sq=lambda v: D(v).sqrt(); ln=A_ln
    f={
     1:lambda: ln(143640)/ln(70), 2:lambda: 3*ln(10)/ln(12), 3:lambda: 3*ln(196)/ln(436), 4:lambda: ln(47216)/ln(48), 5:lambda: ln(140600)/ln(110),
     6:lambda: PI/2, 7:lambda: sq(PI/2), 8:lambda: ln(D(3)/2)/ln(2),
     9:lambda: PI/8, 10:lambda: PI/12, 11:lambda: PI/4,
     12:lambda:(4-sq(6))/5, 13:lambda: sq(4-sq(15)), 14:lambda: 1-1/sq(2),
     15:lambda:(4*PI/11)**D('0.75') if False else (D(3)/4*ln(4*PI/11)).exp(),
     16:lambda: 1+PI*sq(4+1/(PI*PI)), 17:lambda: 2*PI+2,
     18:lambda: D(3)/4-sq(2)/4+sq(6)/2, 19:lambda:(11+sq(61))/12,
     20:lambda: ln((1+sq(5))/2)/2, 21:lambda: ln(1+sq(2))/ln(2), 22:lambda: 3*sq(2*PI)}
    return f[i]()
# ---------- program B ----------
S=10**250
def B_isqrt(x): return isqrt(x*S)          # sqrt of fixed point x
def B_pi():                                  # Gauss-Legendre
    a=S; b=B_isqrt(S//2); t=S//4; p=S
    for _ in range(8):
        an=(a+b)//2; b=B_isqrt(a*b//S); t-=p*(a-an)*(a-an)//S//S; a=an; p*=2
    return (a+b)*(a+b)//S//(4*t)*S//S if False else ((a+b)*(a+b)*S)//(4*t*S//S*1)//1 if False else (a+b)*(a+b)//(4*t)
PIB=B_pi()
def B_agm(a,b):
    for _ in range(60):
        a,b=(a+b)//2,B_isqrt(a*b//S)
    return a
def B_ln(x):                                  # x fixed point; ln via AGM
    m=170
    s=x*(1<<m)                                # s = x*2^m > 2^(p/2)
    ln2=None
    L=PIB*S//(2*B_agm(S,4*S*S//s)) if False else PIB*S//(2*B_agm(S,(4*S*S)//s))
    return L-m*B_LN2
# ln2 by AGM with x=1 (m large)
def _ln2():
    m=240; s=S*(1<<m)*2   # ln(2*2^m)... compute ln(2^m*2) - m ln2 unknown; instead use ln(2^m)=m ln2
    s=S*(1<<m)
    return PIB*S//(2*B_agm(S,(4*S*S)//s))//m
B_LN2=_ln2()
def B_val(i):
    sq=B_isqrt; ln=B_ln
    def e(x):  # exp of fixed point, Taylor
        r=S; t=S; k=0
        while t:
            k+=1; t=t*x//S//k; r+=t
        return r
    LN=lambda n: ln(n*S)
    if i==1: return LN(143640)*S//LN(70)
    if i==2: return 3*LN(10)*S//LN(12)
    if i==3: return 3*LN(196)*S//LN(436)
    if i==4: return LN(47216)*S//LN(48)
    if i==5: return LN(140600)*S//LN(110)
    if i==6: return PIB//2
    if i==7: return sq(PIB//2)
    if i==8: return (LN(3)-LN(2))*S//LN(2)
    if i==9: return PIB//8
    if i==10: return PIB//12
    if i==11: return PIB//4
    if i==12: return (4*S-sq(6*S))//5
    if i==13: return sq(4*S-sq(15*S))
    if i==14: return S-S*S//sq(2*S)
    if i==15: return e(3*(LN(4)+PIB*0)*0 + 3*(ln(4*PIB*S//(11*S)*1))//4) if False else e(3*ln(4*PIB//11)//4)
    if i==16: return S+PIB*sq(4*S+S*S//(PIB*PIB//S))//S
    if i==17: return 2*PIB+2*S
    if i==18: return 3*S//4-sq(2*S)//4+sq(6*S)//2
    if i==19: return (11*S+sq(61*S))//12
    if i==20: return ln((S+sq(5*S))//2)//2
    if i==21: return ln(S+sq(2*S))*S//LN(2)
    if i==22: return 3*sq(2*PIB)
ROWS=[ # id, page, ref, expression, printed, digits(decimals)
 (1,'15a','P1978','log_70 143640','2.7952',4),(2,'15a','BCRL1979,B1980','3 log_12 10','2.7799',4),
 (3,'15a','P1979','3 log_436 196','2.6054',4),(4,'15a','P1980','log_48 47216','2.7802',4),
 (5,'15a','S1981','log_110 140600','2.5218006',7),(6,'11a','BELP2008','pi/2','1.57080',5),
 (7,'11a','Pisier1986','sqrt(pi/2)','1.25331',5),(8,'11b','KP2020','log_2(3/2)','0.58496',5),
 (9,'13a','Wet1973','pi/8','0.3926990816',10),(10,'13a','PW2021','pi/12','0.2617993878',10),
 (11,'13b','trivial','pi/4','0.7853981634',10),(12,'1b','S1958','(4-sqrt6)/5','0.310102',6),
 (13,'1b','M1959','sqrt(4-sqrt15)','0.356393',6),(14,'1b','Scherk1955','1-1/sqrt2','0.292893',6),
 (15,'22a','BS1999','(4pi/11)^(3/4)','1.105',3),(16,'22b','Klotz2021','1+pi sqrt(4+1/pi^2)','7.36',2),
 (17,'22b','trivial','2pi+2','8.28',2),(18,'47a','Ald2000','3/4-sqrt2/4+sqrt6/2','1.6211915',7),
 (19,'47a','Mel2003','(11+sqrt61)/12','1.5675208',7),(20,'82a','Zag93','(1/2)log((1+sqrt5)/2)','0.2406059',7),
 (21,'35a','AP24','log_2(1+sqrt2)','1.271',3),(22,'24a','GFL2026a','3 sqrt(2 pi)','7.5199',4)]
KIND={1: 'u', 2: 'u', 3: 'u', 4: 'u', 5: 'u', 6: 'u', 7: 'l', 8: 'u', 9: 'u', 10: 'u', 11: 'l', 12: 'l', 13: 'l', 14: 'l', 15: 'l', 16: 'u', 17: 'u', 18: 'l', 19: 'l', 20: '?', 21: 'l', 22: 'u'}
def main():
    bad=0
    for (i,pg,ref,ex,pr,nd) in ROWS:
        a=A_val(i); b=D(B_val(i))/D(10**250)
        agree=abs(a-b)<D(10)**-40
        u=D(10)**-nd; p=D(pr)
        ok=abs(a-p)<=u
        red=abs(a-(p+2*u))<=u   # must be False
        near=a.quantize(u)==p
        k=KIND[i]; safe=('n/a' if k=='?' else (p>=a if k=='u' else p<=a))
        st='PASS' if (agree and ok and not red) else 'FAIL'
        if st=='FAIL': bad+=1
        print(f"{st} {pg} [{ref}] {ex} printed {pr} | A={str(a)[:22]} B={str(b)[:22]} agree40={agree} within1ulp={ok} rounds_to_printed={near} red_control_fails={not red} kind={k} printed_on_safe_side={safe}")
    print('ALL PASS' if not bad else f'{bad} FAIL'); sys.exit(1 if bad else 0)
main()
