import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
import sys
from find_core import leader_motif, cut_pep

in_file = sys.argv[1]
out_file = sys.argv[2]



def pep_param(pep):
    
    lanA_param = ProteinAnalysis(pep)
    
    lanA_mw = lanA_param.molecular_weight()
    
    params = [lanA_mw]
    if len(pep) > 0:
        lanA_pI = lanA_param.isoelectric_point()
    else:
        lanA_pI = 'na'   
    params.extend([lanA_pI])
    lanA_AAs = lanA_param.count_amino_acids().values()
    params.extend(lanA_AAs)
    return params

def calc_ratios(core, params):
    S_T = core.count("S") + core.count("T")
    S_T_C = core.count("C") + S_T
    frac_C = float(core.count("C"))/float(len(core))
    frac_S_T = float(S_T)/float(len(core))
    frac_S_T_C = float(S_T_C)/float(len(core))
    charged = core.count("D") + core.count("E") + core.count("K") + core.count("R")
    pos = core.count("K") + core.count("R")
    neg = core.count("D") + core.count("E")
    net_charge = pos - neg
    frac_charged = float(charged)/float(len(core))
    frac_pos = float(pos)/float(len(core))
    frac_neg = float(neg)/float(len(core))
    polar = core.count("C") + core.count("S") + core.count("T") + core.count("H") + core.count("N") + core.count("Q")
    alaphatic = core.count("A") + core.count("I") + core.count("M") + core.count("L") + core.count("V")
    aromatic = core.count("F") + core.count("W") + core.count("Y")
    frac_polar = float(polar)/float(len(core))
    frac_alaphatic = float(alaphatic)/float(len(core))
    frac_aromatic = float(aromatic)/float(len(core))
    params = params + [S_T, S_T_C, frac_C, frac_S_T, frac_S_T_C, charged, pos, neg, net_charge, frac_charged, frac_pos, frac_neg, polar, alaphatic, aromatic, frac_polar, frac_alaphatic, frac_aromatic]
    return params

def positional_fracs(pep, params):
    chunks = int(round(len(pep)/5, 0))
    fracs = []
    for x in range(0, 4):
        fracs.append(float(pep[x*chunks:(x+1)*chunks].count("C"))/float(chunks))
        fracs.append(float((pep[x*chunks:(x+1)*chunks].count("S") + pep[x*chunks:(x+1)*chunks].count("T")))/float(chunks))
        fracs.append(float((pep[x*chunks:(x+1)*chunks].count("S") + pep[x*chunks:(x+1)*chunks].count("T") + pep[x*chunks:(x+1)*chunks].count("C")))/float(chunks))
    fracs.append(float(pep[4*chunks:].count("C"))/float(len(pep[4*chunks:])))
    fracs.append(float((pep[4*chunks:].count("S") + pep[4*chunks:].count("T")))/float(len(pep[4*chunks:])))
    fracs.append(float((pep[4*chunks:].count("S") + pep[4*chunks:].count("T") + pep[4*chunks:].count("C")))/float(len(pep[4*chunks:])))
    params = params + fracs
    return params

def aa_pairs(pep, params, n):
    aa_pairs = []
    aapair_hold = []
    for n in range(1,7):
        pairs = {'AA': 0,'AC': 0,'AE': 0,'AD': 0,'AG': 0,'AF': 0,'AI': 0,'AH': 0,'AK': 0,'AM': 0,'AL': 0,'AN': 0,'AQ': 0,'AP': 0,'AS': 0,'AR': 0,'AT': 0,'AW': 0,'AV': 0,'AY': 0,
'CA': 0,'CC': 0,'CE': 0,'CD': 0,'CG': 0,'CF': 0,'CI': 0,'CH': 0,'CK': 0,'CM': 0,'CL': 0,'CN': 0,'CQ': 0,'CP': 0,'CS': 0,'CR': 0,'CT': 0,'CW': 0,'CV': 0,'CY': 0,
'EA': 0,'EC': 0,'EE': 0,'ED': 0,'EG': 0,'EF': 0,'EI': 0,'EH': 0,'EK': 0,'EM': 0,'EL': 0,'EN': 0,'EQ': 0,'EP': 0,'ES': 0,'ER': 0,'ET': 0,'EW': 0,'EV': 0,'EY': 0,
'DA': 0,'DC': 0,'DE': 0,'DD': 0,'DG': 0,'DF': 0,'DI': 0,'DH': 0,'DK': 0,'DM': 0,'DL': 0,'DN': 0,'DQ': 0,'DP': 0,'DS': 0,'DR': 0,'DT': 0,'DW': 0,'DV': 0,'DY': 0,
'GA': 0,'GC': 0,'GE': 0,'GD': 0,'GG': 0,'GF': 0,'GI': 0,'GH': 0,'GK': 0,'GM': 0,'GL': 0,'GN': 0,'GQ': 0,'GP': 0,'GS': 0,'GR': 0,'GT': 0,'GW': 0,'GV': 0,'GY': 0,
'FA': 0,'FC': 0,'FE': 0,'FD': 0,'FG': 0,'FF': 0,'FI': 0,'FH': 0,'FK': 0,'FM': 0,'FL': 0,'FN': 0,'FQ': 0,'FP': 0,'FS': 0,'FR': 0,'FT': 0,'FW': 0,'FV': 0,'FY': 0,
'IA': 0,'IC': 0,'IE': 0,'ID': 0,'IG': 0,'IF': 0,'II': 0,'IH': 0,'IK': 0,'IM': 0,'IL': 0,'IN': 0,'IQ': 0,'IP': 0,'IS': 0,'IR': 0,'IT': 0,'IW': 0,'IV': 0,'IY': 0,
'HA': 0,'HC': 0,'HE': 0,'HD': 0,'HG': 0,'HF': 0,'HI': 0,'HH': 0,'HK': 0,'HM': 0,'HL': 0,'HN': 0,'HQ': 0,'HP': 0,'HS': 0,'HR': 0,'HT': 0,'HW': 0,'HV': 0,'HY': 0,
'KA': 0,'KC': 0,'KE': 0,'KD': 0,'KG': 0,'KF': 0,'KI': 0,'KH': 0,'KK': 0,'KM': 0,'KL': 0,'KN': 0,'KQ': 0,'KP': 0,'KS': 0,'KR': 0,'KT': 0,'KW': 0,'KV': 0,'KY': 0,
'MA': 0,'MC': 0,'ME': 0,'MD': 0,'MG': 0,'MF': 0,'MI': 0,'MH': 0,'MK': 0,'MM': 0,'ML': 0,'MN': 0,'MQ': 0,'MP': 0,'MS': 0,'MR': 0,'MT': 0,'MW': 0,'MV': 0,'MY': 0,
'LA': 0,'LC': 0,'LE': 0,'LD': 0,'LG': 0,'LF': 0,'LI': 0,'LH': 0,'LK': 0,'LM': 0,'LL': 0,'LN': 0,'LQ': 0,'LP': 0,'LS': 0,'LR': 0,'LT': 0,'LW': 0,'LV': 0,'LY': 0,
'NA': 0,'NC': 0,'NE': 0,'ND': 0,'NG': 0,'NF': 0,'NI': 0,'NH': 0,'NK': 0,'NM': 0,'NL': 0,'NN': 0,'NQ': 0,'NP': 0,'NS': 0,'NR': 0,'NT': 0,'NW': 0,'NV': 0,'NY': 0,
'QA': 0,'QC': 0,'QE': 0,'QD': 0,'QG': 0,'QF': 0,'QI': 0,'QH': 0,'QK': 0,'QM': 0,'QL': 0,'QN': 0,'QQ': 0,'QP': 0,'QS': 0,'QR': 0,'QT': 0,'QW': 0,'QV': 0,'QY': 0,
'PA': 0,'PC': 0,'PE': 0,'PD': 0,'PG': 0,'PF': 0,'PI': 0,'PH': 0,'PK': 0,'PM': 0,'PL': 0,'PN': 0,'PQ': 0,'PP': 0,'PS': 0,'PR': 0,'PT': 0,'PW': 0,'PV': 0,'PY': 0,
'SA': 0,'SC': 0,'SE': 0,'SD': 0,'SG': 0,'SF': 0,'SI': 0,'SH': 0,'SK': 0,'SM': 0,'SL': 0,'SN': 0,'SQ': 0,'SP': 0,'SS': 0,'SR': 0,'ST': 0,'SW': 0,'SV': 0,'SY': 0,
'RA': 0,'RC': 0,'RE': 0,'RD': 0,'RG': 0,'RF': 0,'RI': 0,'RH': 0,'RK': 0,'RM': 0,'RL': 0,'RN': 0,'RQ': 0,'RP': 0,'RS': 0,'RR': 0,'RT': 0,'RW': 0,'RV': 0,'RY': 0,
'TA': 0,'TC': 0,'TE': 0,'TD': 0,'TG': 0,'TF': 0,'TI': 0,'TH': 0,'TK': 0,'TM': 0,'TL': 0,'TN': 0,'TQ': 0,'TP': 0,'TS': 0,'TR': 0,'TT': 0,'TW': 0,'TV': 0,'TY': 0,
'WA': 0,'WC': 0,'WE': 0,'WD': 0,'WG': 0,'WF': 0,'WI': 0,'WH': 0,'WK': 0,'WM': 0,'WL': 0,'WN': 0,'WQ': 0,'WP': 0,'WS': 0,'WR': 0,'WT': 0,'WW': 0,'WV': 0,'WY': 0,
'VA': 0,'VC': 0,'VE': 0,'VD': 0,'VG': 0,'VF': 0,'VI': 0,'VH': 0,'VK': 0,'VM': 0,'VL': 0,'VN': 0,'VQ': 0,'VP': 0,'VS': 0,'VR': 0,'VT': 0,'VW': 0,'VV': 0,'VY': 0,
'YA': 0,'YC': 0,'YE': 0,'YD': 0,'YG': 0,'YF': 0,'YI': 0,'YH': 0,'YK': 0,'YM': 0,'YL': 0,'YN': 0,'YQ': 0,'YP': 0,'YS': 0,'YR': 0,'YT': 0,'YW': 0,'YV': 0,'YY': 0}
    
        for x in range(0, len(pep)-n):
            pairs[pep[x]+pep[x+n]] = pairs[pep[x]+pep[x+n]] + 1

        for res1 in ['A','C','E','D','G','F','I','H','K','M','L','N','Q','P','S','R','T','W','V','Y']:
            for res2 in ['A','C','E','D','G','F','I','H','K','M','L','N','Q','P','S','R','T','W','V','Y']:
                aapair_hold.append(pairs[res1+res2])
    for i in [1,4,7,10,9,13,15,16,18,17,400,401,407,410,411,413,415,416,418,417,800,804,807,810,813,815,816,1200,1207,1210,1211,1213,1215,1600,1601,1610,1609,1613,1615,
              1614,1617,2000,2010,2011,2013,2015,2016,2017,20,21,23,22,24,26,28,31,33,35,34,36,38,37,421,423,422,424,428,431,433,435,434,436,438,821,823,824,828,831,832,
              835,834,836,838,1221,1223,1224,1227,1226,1228,1231,1235,1234,1236,1620,1621,1623,1624,1626,1631,1635,1634,1636,2021,2025,2024,2026,2028,2031,2032,2035,2034,
              2036,839,61,63,65,64,67,71,73,75,74,76,461,462,464,471,475,476,861,864,871,874,876,1261,1263,1264,1268,1271,1275,1274,1276,1661,1670,1671,1675,1674,1676,2061,
              2063,2062,2071,2075,2076,1679,43,55,56,455,456,841,855,856,1243,1251,1655,1656,2041,2047,2051,2055,2056,100,105,106,110,113,112,115,116,500,501,505,506,510,513,
              515,901,910,911,915,914,1301,1305,1306,1310,1315,1314,1318,1317,1701,1704,1706,1710,1715,2105,2106,2108,2110,2115,2116,81,83,82,87,90,91,93,95,94,96,97,480,481,
              490,491,493,495,494,496,498,880,881,886,890,891,893,895,896,898,1280,1281,1288,1290,1291,1293,1295,1294,1296,1680,1681,1687,1688,1690,1691,1693,1695,1694,1696,
              2081,2083,2090,2091,2093,2095,2094,2096,140,141,150,153,155,540,544,555,940,941,944,950,954,1340,1341,1350,1353,1355,1354,1740,1750,1755,2148,2150,2153,2155,121,
              125,135,136,521,530,535,534,536,920,921,925,936,1335,1334,1336,1721,1725,1735,1736,2121,2125,2130,2135,2136,161,176,561,564,570,575,574,576,960,961,963,975,976,
              1361,1362,1370,1376,1761,1763,1762,1764,1775,1776,2161,2164,2170,2175,2174,200,205,207,208,210,213,215,214,218,605,604,607,610,609,611,613,612,615,614,618,1005,
              1004,1007,1008,1010,1013,1015,1014,1018,1400,1405,1407,1410,1413,1415,1414,1800,1801,1804,1807,1806,1810,1813,1815,1814,1818,2200,2205,2207,2206,2210,2209,2212,
              2215,2218,2217,2219,619,181,198,197,595,993,995,996,1395,1795,221,224,235,234,236,621,623,624,628,631,635,634,636,1021,1024,1032,1035,1034,1036,1038,1421,1424,
              1433,1436,1820,1821,1822,1831,1833,1832,1836,2221,2223,2224,2234,2236,260,263,265,267,270,271,273,275,274,276,278,277,660,661,664,670,673,675,674,676,678,677,
              1060,1070,1073,1075,1074,1460,1464,1470,1473,1475,1474,1476,1860,1864,1870,1873,1875,1878,1877,2260,2261,2263,2264,2270,2273,2275,2278,241,255,256,655,1041,1045,
              1055,1450,1455,1841,1850,1855,2241,2247,2250,2255,300,301,303,302,305,304,307,306,308,310,309,311,313,312,315,314,316,318,317,700,701,703,702,705,704,707,708,710,
              711,713,712,715,714,716,718,717,1100,1103,1105,1104,1107,1106,1110,1109,1113,1112,1115,1114,1116,1118,1117,1500,1501,1503,1502,1505,1504,1507,1506,1510,1509,1513,
              1512,1515,1514,1516,1518,1517,1900,1901,1903,1902,1905,1904,1907,1908,1910,1909,1913,1912,1915,1914,1916,1918,1917,2300,2301,2303,2302,2305,2304,2307,2306,2310,
              2309,2311,2313,2312,2315,2314,2316,2318,2317,2319,1119,719,280,281,283,284,288,291,293,295,294,296,680,681,683,684,688,690,689,691,693,695,694,697,1081,1083,1088,
              1090,1093,1095,1096,1480,1481,1483,1490,1491,1495,1496,1497,1880,1881,1888,1890,1891,1893,1895,1896,1897,2280,2281,2290,2291,2295,2296,2297,321,323,322,324,327,326,
              328,331,332,334,336,720,721,723,722,724,726,728,731,735,734,736,738,1120,1121,1123,1125,1124,1126,1128,1131,1135,1134,1136,1138,1520,1521,1523,1522,1525,1524,1526,
              1528,1531,1532,1535,1534,1536,1921,1923,1922,1924,1926,1928,1931,1932,1935,1934,1936,2321,2322,2324,2328,2331,2335,2334,2336,739,361,365,364,366,370,373,375,376,377,
              760,761,764,768,770,773,775,778,1164,1167,1170,1175,1176,1178,1560,1564,1570,1573,1575,1970,1973,1975,1978,2360,2361,2370,2373,2375,2374,2378,340,353,355,740,753,
              755,1153,1155,1154,1540,1553,1555,1940,1955,2340,2350,2355,381,390,395,396,781,795,1181,1190,1590,1996,2395]:
        aa_pairs.append(aapair_hold[i])
    return params




fh = open(in_file)
out = open(out_file, 'w')

csv.writer(out).writerow(['query_acc', 'phylum', 'genus/species', 'nucleotide_acc', 'min', 'max', 'dir', 'peptide', 'predicted_leader', 'predicted_core', 'peptide_length', 
    'leader_length', 'core_length', 'core_MW','core_pI','core_A','core_C','core_E','core_D','core_G','core_F','core_I','core_H','core_K','core_M','core_L','core_N','core_Q',
    'core_P','core_S','core_R','core_T','core_W','core_V','core_Y',
    "S+T", "S+T+C", "frac_c", "frac_S+T", "frac_S+T+C",  "charged", "pos", "neg", "net_charge", "frac_charge", "frac_pos", "frac_neg", "polar", "alaphatic", "aromatic", "frac_polar", "frac_alaphatic", "frac_aromatic",
    "1/5_frac_c", "1/5_frac_s_t", "1/5_frac_s_t_c","2/5_frac_c", "2/5_frac_s_t", "2/5_frac_s_t_c","3/5_frac_c", "3/5_frac_s_t", "3/5_frac_s_t_c","4/5_frac_c", "4/5_frac_s_t", "4/5_frac_s_t_c","5/5_frac_c", "5/5_frac_s_t", "5/5_frac_s_t_c",
    'core_AA','core_AC','core_AE','core_AD','core_AG','core_AF','core_AI','core_AH','core_AK','core_AM','core_AL',
    'core_AN','core_AQ','core_AP','core_AS','core_AR','core_AT','core_AW','core_AV','core_AY','core_CA','core_CC','core_CE','core_CD','core_CG','core_CF','core_CI','core_CH',
    'core_CK','core_CM','core_CL','core_CN','core_CQ','core_CP','core_CS','core_CR','core_CT','core_CW','core_CV','core_CY','core_EA','core_EC','core_EE','core_ED','core_EG',
    'core_EF','core_EI','core_EH','core_EK','core_EM','core_EL','core_EN','core_EQ','core_EP','core_ES','core_ER','core_ET','core_EW','core_EV','core_EY','core_DA','core_DC',
    'core_DE','core_DD','core_DG','core_DF','core_DI','core_DH','core_DK','core_DM','core_DL','core_DN','core_DQ','core_DP','core_DS','core_DR','core_DT','core_DW','core_DV',
    'core_DY','core_GA','core_GC','core_GE','core_GD','core_GG','core_GF','core_GI','core_GH','core_GK','core_GM','core_GL','core_GN','core_GQ','core_GP','core_GS','core_GR',
    'core_GT','core_GW','core_GV','core_GY','core_FA','core_FC','core_FE','core_FD','core_FG','core_FF','core_FI','core_FH','core_FK','core_FM','core_FL','core_FN','core_FQ',
    'core_FP','core_FS','core_FR','core_FT','core_FW','core_FV','core_FY','core_IA','core_IC','core_IE','core_ID','core_IG','core_IF','core_II','core_IH','core_IK','core_IM',
    'core_IL','core_IN','core_IQ','core_IP','core_IS','core_IR','core_IT','core_IW','core_IV','core_IY','core_HA','core_HC','core_HE','core_HD','core_HG','core_HF','core_HI',
    'core_HH','core_HK','core_HM','core_HL','core_HN','core_HQ','core_HP','core_HS','core_HR','core_HT','core_HW','core_HV','core_HY','core_KA','core_KC','core_KE','core_KD',
    'core_KG','core_KF','core_KI','core_KH','core_KK','core_KM','core_KL','core_KN','core_KQ','core_KP','core_KS','core_KR','core_KT','core_KW','core_KV','core_KY','core_MA',
    'core_MC','core_ME','core_MD','core_MG','core_MF','core_MI','core_MH','core_MK','core_MM','core_ML','core_MN','core_MQ','core_MP','core_MS','core_MR','core_MT','core_MW',
    'core_MV','core_MY','core_LA','core_LC','core_LE','core_LD','core_LG','core_LF','core_LI','core_LH','core_LK','core_LM','core_LL','core_LN','core_LQ','core_LP','core_LS',
    'core_LR','core_LT','core_LW','core_LV','core_LY','core_NA','core_NC','core_NE','core_ND','core_NG','core_NF','core_NI','core_NH','core_NK','core_NM','core_NL','core_NN',
    'core_NQ','core_NP','core_NS','core_NR','core_NT','core_NW','core_NV','core_NY','core_QA','core_QC','core_QE','core_QD','core_QG','core_QF','core_QI','core_QH','core_QK',
    'core_QM','core_QL','core_QN','core_QQ','core_QP','core_QS','core_QR','core_QT','core_QW','core_QV','core_QY','core_PA','core_PC','core_PE','core_PD','core_PG','core_PF',
    'core_PI','core_PH','core_PK','core_PM','core_PL','core_PN','core_PQ','core_PP','core_PS','core_PR','core_PT','core_PW','core_PV','core_PY','core_SA','core_SC','core_SE',
    'core_SD','core_SG','core_SF','core_SI','core_SH','core_SK','core_SM','core_SL','core_SN','core_SQ','core_SP','core_SS','core_SR','core_ST','core_SW','core_SV','core_SY',
    'core_RA','core_RC','core_RE','core_RD','core_RG','core_RF','core_RI','core_RH','core_RK','core_RM','core_RL','core_RN','core_RQ','core_RP','core_RS','core_RR','core_RT',
    'core_RW','core_RV','core_RY','core_TA','core_TC','core_TE','core_TD','core_TG','core_TF','core_TI','core_TH','core_TK','core_TM','core_TL','core_TN','core_TQ','core_TP',
    'core_TS','core_TR','core_TT','core_TW','core_TV','core_TY','core_WA','core_WC','core_WE','core_WD','core_WG','core_WF','core_WI','core_WH','core_WK','core_WM','core_WL',
    'core_WN','core_WQ','core_WP','core_WS','core_WR','core_WT','core_WW','core_WV','core_WY','core_VA','core_VC','core_VE','core_VD','core_VG','core_VF','core_VI','core_VH',
    'core_VK','core_VM','core_VL','core_VN','core_VQ','core_VP','core_VS','core_VR','core_VT','core_VW','core_VV','core_VY','core_YA','core_YC','core_YE','core_YD','core_YG',
    'core_YF','core_YI','core_YH','core_YK','core_YM','core_YL','core_YN','core_YQ','core_YP','core_YS','core_YR','core_YT','core_YW','core_YV','core_YY','core_AxA','core_AxC',
    'core_AxE','core_AxD','core_AxG','core_AxF','core_AxI','core_AxH','core_AxK','core_AxM','core_AxL','core_AxN','core_AxQ','core_AxP','core_AxS','core_AxR','core_AxT','core_AxW',
    'core_AxV','core_AxY','core_CxA','core_CxC','core_CxE','core_CxD','core_CxG','core_CxF','core_CxI','core_CxH','core_CxK','core_CxM','core_CxL','core_CxN','core_CxQ','core_CxP',
    'core_CxS','core_CxR','core_CxT','core_CxW','core_CxV','core_CxY','core_ExA','core_ExC','core_ExE','core_ExD','core_ExG','core_ExF','core_ExI','core_ExH','core_ExK','core_ExM',
    'core_ExL','core_ExN','core_ExQ','core_ExP','core_ExS','core_ExR','core_ExT','core_ExW','core_ExV','core_ExY','core_DxA','core_DxC','core_DxE','core_DxD','core_DxG','core_DxF',
    'core_DxI','core_DxH','core_DxK','core_DxM','core_DxL','core_DxN','core_DxQ','core_DxP','core_DxS','core_DxR','core_DxT','core_DxW','core_DxV','core_DxY','core_GxA','core_GxC',
    'core_GxE','core_GxD','core_GxG','core_GxF','core_GxI','core_GxH','core_GxK','core_GxM','core_GxL','core_GxN','core_GxQ','core_GxP','core_GxS','core_GxR','core_GxT','core_GxW',
    'core_GxV','core_GxY','core_FxA','core_FxC','core_FxE','core_FxD','core_FxG','core_FxF','core_FxI','core_FxH','core_FxK','core_FxM','core_FxL','core_FxN','core_FxQ','core_FxP',
    'core_FxS','core_FxR','core_FxT','core_FxW','core_FxV','core_FxY','core_IxA','core_IxC','core_IxE','core_IxD','core_IxG','core_IxF','core_IxI','core_IxH','core_IxK','core_IxM',
    'core_IxL','core_IxN','core_IxQ','core_IxP','core_IxS','core_IxR','core_IxT','core_IxW','core_IxV','core_IxY','core_HxA','core_HxC','core_HxE','core_HxD','core_HxG','core_HxF',
    'core_HxI','core_HxH','core_HxK','core_HxM','core_HxL','core_HxN','core_HxQ','core_HxP','core_HxS','core_HxR','core_HxT','core_HxW','core_HxV','core_HxY','core_KxA','core_KxC',
    'core_KxE','core_KxD','core_KxG','core_KxF','core_KxI','core_KxH','core_KxK','core_KxM','core_KxL','core_KxN','core_KxQ','core_KxP','core_KxS','core_KxR','core_KxT','core_KxW',
    'core_KxV','core_KxY','core_MxA','core_MxC','core_MxE','core_MxD','core_MxG','core_MxF','core_MxI','core_MxH','core_MxK','core_MxM','core_MxL','core_MxN','core_MxQ','core_MxP',
    'core_MxS','core_MxR','core_MxT','core_MxW','core_MxV','core_MxY','core_LxA','core_LxC','core_LxE','core_LxD','core_LxG','core_LxF','core_LxI','core_LxH','core_LxK','core_LxM',
    'core_LxL','core_LxN','core_LxQ','core_LxP','core_LxS','core_LxR','core_LxT','core_LxW','core_LxV','core_LxY','core_NxA','core_NxC','core_NxE','core_NxD','core_NxG','core_NxF',
    'core_NxI','core_NxH','core_NxK','core_NxM','core_NxL','core_NxN','core_NxQ','core_NxP','core_NxS','core_NxR','core_NxT','core_NxW','core_NxV','core_NxY','core_QxA','core_QxC',
    'core_QxE','core_QxD','core_QxG','core_QxF','core_QxI','core_QxH','core_QxK','core_QxM','core_QxL','core_QxN','core_QxQ','core_QxP','core_QxS','core_QxR','core_QxT','core_QxW',
    'core_QxV','core_QxY','core_PxA','core_PxC','core_PxE','core_PxD','core_PxG','core_PxF','core_PxI','core_PxH','core_PxK','core_PxM','core_PxL','core_PxN','core_PxQ','core_PxP',
    'core_PxS','core_PxR','core_PxT','core_PxW','core_PxV','core_PxY','core_SxA','core_SxC','core_SxE','core_SxD','core_SxG','core_SxF','core_SxI','core_SxH','core_SxK','core_SxM',
    'core_SxL','core_SxN','core_SxQ','core_SxP','core_SxS','core_SxR','core_SxT','core_SxW','core_SxV','core_SxY','core_RxA','core_RxC','core_RxE','core_RxD','core_RxG','core_RxF',
    'core_RxI','core_RxH','core_RxK','core_RxM','core_RxL','core_RxN','core_RxQ','core_RxP','core_RxS','core_RxR','core_RxT','core_RxW','core_RxV','core_RxY','core_TxA','core_TxC',
    'core_TxE','core_TxD','core_TxG','core_TxF','core_TxI','core_TxH','core_TxK','core_TxM','core_TxL','core_TxN','core_TxQ','core_TxP','core_TxS','core_TxR','core_TxT','core_TxW',
    'core_TxV','core_TxY','core_WxA','core_WxC','core_WxE','core_WxD','core_WxG','core_WxF','core_WxI','core_WxH','core_WxK','core_WxM','core_WxL','core_WxN','core_WxQ','core_WxP',
    'core_WxS','core_WxR','core_WxT','core_WxW','core_WxV','core_WxY','core_VxA','core_VxC','core_VxE','core_VxD','core_VxG','core_VxF','core_VxI','core_VxH','core_VxK','core_VxM',
    'core_VxL','core_VxN','core_VxQ','core_VxP','core_VxS','core_VxR','core_VxT','core_VxW','core_VxV','core_VxY','core_YxA','core_YxC','core_YxE','core_YxD','core_YxG','core_YxF',
    'core_YxI','core_YxH','core_YxK','core_YxM','core_YxL','core_YxN','core_YxQ','core_YxP','core_YxS','core_YxR','core_YxT','core_YxW','core_YxV','core_YxY','core_AxxA','core_AxxC',
'core_AxxE','core_AxxD','core_AxxG','core_AxxF','core_AxxI','core_AxxH','core_AxxK','core_AxxM','core_AxxL','core_AxxN','core_AxxQ','core_AxxP','core_AxxS','core_AxxR','core_AxxT',
'core_AxxW','core_AxxV','core_AxxY','core_CxxA','core_CxxC','core_CxxE','core_CxxD','core_CxxG','core_CxxF','core_CxxI','core_CxxH','core_CxxK','core_CxxM','core_CxxL','core_CxxN',
'core_CxxQ','core_CxxP','core_CxxS','core_CxxR','core_CxxT','core_CxxW','core_CxxV','core_CxxY','core_ExxA','core_ExxC','core_ExxE','core_ExxD','core_ExxG','core_ExxF','core_ExxI',
'core_ExxH','core_ExxK','core_ExxM','core_ExxL','core_ExxN','core_ExxQ','core_ExxP','core_ExxS','core_ExxR','core_ExxT','core_ExxW','core_ExxV','core_ExxY','core_DxxA','core_DxxC',
'core_DxxE','core_DxxD','core_DxxG','core_DxxF','core_DxxI','core_DxxH','core_DxxK','core_DxxM','core_DxxL','core_DxxN','core_DxxQ','core_DxxP','core_DxxS','core_DxxR','core_DxxT',
'core_DxxW','core_DxxV','core_DxxY','core_GxxA','core_GxxC','core_GxxE','core_GxxD','core_GxxG','core_GxxF','core_GxxI','core_GxxH','core_GxxK','core_GxxM','core_GxxL','core_GxxN',
'core_GxxQ','core_GxxP','core_GxxS','core_GxxR','core_GxxT','core_GxxW','core_GxxV','core_GxxY','core_FxxA','core_FxxC','core_FxxE','core_FxxD','core_FxxG','core_FxxF','core_FxxI',
'core_FxxH','core_FxxK','core_FxxM','core_FxxL','core_FxxN','core_FxxQ','core_FxxP','core_FxxS','core_FxxR','core_FxxT','core_FxxW','core_FxxV','core_FxxY','core_IxxA','core_IxxC',
'core_IxxE','core_IxxD','core_IxxG','core_IxxF','core_IxxI','core_IxxH','core_IxxK','core_IxxM','core_IxxL','core_IxxN','core_IxxQ','core_IxxP','core_IxxS','core_IxxR','core_IxxT',
'core_IxxW','core_IxxV','core_IxxY','core_HxxA','core_HxxC','core_HxxE','core_HxxD','core_HxxG','core_HxxF','core_HxxI','core_HxxH','core_HxxK','core_HxxM','core_HxxL','core_HxxN',
'core_HxxQ','core_HxxP','core_HxxS','core_HxxR','core_HxxT','core_HxxW','core_HxxV','core_HxxY','core_KxxA','core_KxxC','core_KxxE','core_KxxD','core_KxxG','core_KxxF','core_KxxI',
'core_KxxH','core_KxxK','core_KxxM','core_KxxL','core_KxxN','core_KxxQ','core_KxxP','core_KxxS','core_KxxR','core_KxxT','core_KxxW','core_KxxV','core_KxxY','core_MxxA','core_MxxC',
'core_MxxE','core_MxxD','core_MxxG','core_MxxF','core_MxxI','core_MxxH','core_MxxK','core_MxxM','core_MxxL','core_MxxN','core_MxxQ','core_MxxP','core_MxxS','core_MxxR','core_MxxT',
'core_MxxW','core_MxxV','core_MxxY','core_LxxA','core_LxxC','core_LxxE','core_LxxD','core_LxxG','core_LxxF','core_LxxI','core_LxxH','core_LxxK','core_LxxM','core_LxxL','core_LxxN',
'core_LxxQ','core_LxxP','core_LxxS','core_LxxR','core_LxxT','core_LxxW','core_LxxV','core_LxxY','core_NxxA','core_NxxC','core_NxxE','core_NxxD','core_NxxG','core_NxxF','core_NxxI',
'core_NxxH','core_NxxK','core_NxxM','core_NxxL','core_NxxN','core_NxxQ','core_NxxP','core_NxxS','core_NxxR','core_NxxT','core_NxxW','core_NxxV','core_NxxY','core_QxxA','core_QxxC',
'core_QxxE','core_QxxD','core_QxxG','core_QxxF','core_QxxI','core_QxxH','core_QxxK','core_QxxM','core_QxxL','core_QxxN','core_QxxQ','core_QxxP','core_QxxS','core_QxxR','core_QxxT',
'core_QxxW','core_QxxV','core_QxxY','core_PxxA','core_PxxC','core_PxxE','core_PxxD','core_PxxG','core_PxxF','core_PxxI','core_PxxH','core_PxxK','core_PxxM','core_PxxL','core_PxxN',
'core_PxxQ','core_PxxP','core_PxxS','core_PxxR','core_PxxT','core_PxxW','core_PxxV','core_PxxY','core_SxxA','core_SxxC','core_SxxE','core_SxxD','core_SxxG','core_SxxF','core_SxxI',
'core_SxxH','core_SxxK','core_SxxM','core_SxxL','core_SxxN','core_SxxQ','core_SxxP','core_SxxS','core_SxxR','core_SxxT','core_SxxW','core_SxxV','core_SxxY','core_RxxA','core_RxxC',
'core_RxxE','core_RxxD','core_RxxG','core_RxxF','core_RxxI','core_RxxH','core_RxxK','core_RxxM','core_RxxL','core_RxxN','core_RxxQ','core_RxxP','core_RxxS','core_RxxR','core_RxxT',
'core_RxxW','core_RxxV','core_RxxY','core_TxxA','core_TxxC','core_TxxE','core_TxxD','core_TxxG','core_TxxF','core_TxxI','core_TxxH','core_TxxK','core_TxxM','core_TxxL','core_TxxN',
'core_TxxQ','core_TxxP','core_TxxS','core_TxxR','core_TxxT','core_TxxW','core_TxxV','core_TxxY','core_WxxA','core_WxxC','core_WxxE','core_WxxD','core_WxxG','core_WxxF','core_WxxI',
'core_WxxH','core_WxxK','core_WxxM','core_WxxL','core_WxxN','core_WxxQ','core_WxxP','core_WxxS','core_WxxR','core_WxxT','core_WxxW','core_WxxV','core_WxxY','core_VxxA','core_VxxC',
'core_VxxE','core_VxxD','core_VxxG','core_VxxF','core_VxxI','core_VxxH','core_VxxK','core_VxxM','core_VxxL','core_VxxN','core_VxxQ','core_VxxP','core_VxxS','core_VxxR','core_VxxT',
'core_VxxW','core_VxxV','core_VxxY','core_YxxA','core_YxxC','core_YxxE','core_YxxD','core_YxxG','core_YxxF','core_YxxI','core_YxxH','core_YxxK','core_YxxM','core_YxxL','core_YxxN',
'core_YxxQ','core_YxxP','core_YxxS','core_YxxR','core_YxxT','core_YxxW','core_YxxV','core_YxxY','core_AxxxA','core_AxxxC','core_AxxxE','core_AxxxD','core_AxxxG','core_AxxxF',
'core_AxxxI','core_AxxxH','core_AxxxK','core_AxxxM','core_AxxxL','core_AxxxN','core_AxxxQ','core_AxxxP','core_AxxxS','core_AxxxR','core_AxxxT','core_AxxxW','core_AxxxV','core_AxxxY',
'core_CxxxA','core_CxxxC','core_CxxxE','core_CxxxD','core_CxxxG','core_CxxxF','core_CxxxI','core_CxxxH','core_CxxxK','core_CxxxM','core_CxxxL','core_CxxxN','core_CxxxQ','core_CxxxP',
'core_CxxxS','core_CxxxR','core_CxxxT','core_CxxxW','core_CxxxV','core_CxxxY','core_ExxxA','core_ExxxC','core_ExxxE','core_ExxxD','core_ExxxG','core_ExxxF','core_ExxxI','core_ExxxH',
'core_ExxxK','core_ExxxM','core_ExxxL','core_ExxxN','core_ExxxQ','core_ExxxP','core_ExxxS','core_ExxxR','core_ExxxT','core_ExxxW','core_ExxxV','core_ExxxY','core_DxxxA','core_DxxxC',
'core_DxxxE','core_DxxxD','core_DxxxG','core_DxxxF','core_DxxxI','core_DxxxH','core_DxxxK','core_DxxxM','core_DxxxL','core_DxxxN','core_DxxxQ','core_DxxxP','core_DxxxS','core_DxxxR',
'core_DxxxT','core_DxxxW','core_DxxxV','core_DxxxY','core_GxxxA','core_GxxxC','core_GxxxE','core_GxxxD','core_GxxxG','core_GxxxF','core_GxxxI','core_GxxxH','core_GxxxK','core_GxxxM',
'core_GxxxL','core_GxxxN','core_GxxxQ','core_GxxxP','core_GxxxS','core_GxxxR','core_GxxxT','core_GxxxW','core_GxxxV','core_GxxxY','core_FxxxA','core_FxxxC','core_FxxxE','core_FxxxD',
'core_FxxxG','core_FxxxF','core_FxxxI','core_FxxxH','core_FxxxK','core_FxxxM','core_FxxxL','core_FxxxN','core_FxxxQ','core_FxxxP','core_FxxxS','core_FxxxR','core_FxxxT','core_FxxxW',
'core_FxxxV','core_FxxxY','core_IxxxA','core_IxxxC','core_IxxxE','core_IxxxD','core_IxxxG','core_IxxxF','core_IxxxI','core_IxxxH','core_IxxxK','core_IxxxM','core_IxxxL','core_IxxxN',
'core_IxxxQ','core_IxxxP','core_IxxxS','core_IxxxR','core_IxxxT','core_IxxxW','core_IxxxV','core_IxxxY','core_HxxxA','core_HxxxC','core_HxxxE','core_HxxxD','core_HxxxG','core_HxxxF',
'core_HxxxI','core_HxxxH','core_HxxxK','core_HxxxM','core_HxxxL','core_HxxxN','core_HxxxQ','core_HxxxP','core_HxxxS','core_HxxxR','core_HxxxT','core_HxxxW','core_HxxxV','core_HxxxY',
'core_KxxxA','core_KxxxC','core_KxxxE','core_KxxxD','core_KxxxG','core_KxxxF','core_KxxxI','core_KxxxH','core_KxxxK','core_KxxxM','core_KxxxL','core_KxxxN','core_KxxxQ','core_KxxxP',
'core_KxxxS','core_KxxxR','core_KxxxT','core_KxxxW','core_KxxxV','core_KxxxY','core_MxxxA','core_MxxxC','core_MxxxE','core_MxxxD','core_MxxxG','core_MxxxF','core_MxxxI','core_MxxxH',
'core_MxxxK','core_MxxxM','core_MxxxL','core_MxxxN','core_MxxxQ','core_MxxxP','core_MxxxS','core_MxxxR','core_MxxxT','core_MxxxW','core_MxxxV','core_MxxxY','core_LxxxA','core_LxxxC',
'core_LxxxE','core_LxxxD','core_LxxxG','core_LxxxF','core_LxxxI','core_LxxxH','core_LxxxK','core_LxxxM','core_LxxxL','core_LxxxN','core_LxxxQ','core_LxxxP','core_LxxxS','core_LxxxR',
'core_LxxxT','core_LxxxW','core_LxxxV','core_LxxxY','core_NxxxA','core_NxxxC','core_NxxxE','core_NxxxD','core_NxxxG','core_NxxxF','core_NxxxI','core_NxxxH','core_NxxxK','core_NxxxM',
'core_NxxxL','core_NxxxN','core_NxxxQ','core_NxxxP','core_NxxxS','core_NxxxR','core_NxxxT','core_NxxxW','core_NxxxV','core_NxxxY','core_QxxxA','core_QxxxC','core_QxxxE','core_QxxxD',
'core_QxxxG','core_QxxxF','core_QxxxI','core_QxxxH','core_QxxxK','core_QxxxM','core_QxxxL','core_QxxxN','core_QxxxQ','core_QxxxP','core_QxxxS','core_QxxxR','core_QxxxT','core_QxxxW',
'core_QxxxV','core_QxxxY','core_PxxxA','core_PxxxC','core_PxxxE','core_PxxxD','core_PxxxG','core_PxxxF','core_PxxxI','core_PxxxH','core_PxxxK','core_PxxxM','core_PxxxL','core_PxxxN',
'core_PxxxQ','core_PxxxP','core_PxxxS','core_PxxxR','core_PxxxT','core_PxxxW','core_PxxxV','core_PxxxY','core_SxxxA','core_SxxxC','core_SxxxE','core_SxxxD','core_SxxxG','core_SxxxF',
'core_SxxxI','core_SxxxH','core_SxxxK','core_SxxxM','core_SxxxL','core_SxxxN','core_SxxxQ','core_SxxxP','core_SxxxS','core_SxxxR','core_SxxxT','core_SxxxW','core_SxxxV','core_SxxxY',
'core_RxxxA','core_RxxxC','core_RxxxE','core_RxxxD','core_RxxxG','core_RxxxF','core_RxxxI','core_RxxxH','core_RxxxK','core_RxxxM','core_RxxxL','core_RxxxN','core_RxxxQ','core_RxxxP',
'core_RxxxS','core_RxxxR','core_RxxxT','core_RxxxW','core_RxxxV','core_RxxxY','core_TxxxA','core_TxxxC','core_TxxxE','core_TxxxD','core_TxxxG','core_TxxxF','core_TxxxI','core_TxxxH',
'core_TxxxK','core_TxxxM','core_TxxxL','core_TxxxN','core_TxxxQ','core_TxxxP','core_TxxxS','core_TxxxR','core_TxxxT','core_TxxxW','core_TxxxV','core_TxxxY','core_WxxxA','core_WxxxC',
'core_WxxxE','core_WxxxD','core_WxxxG','core_WxxxF','core_WxxxI','core_WxxxH','core_WxxxK','core_WxxxM','core_WxxxL','core_WxxxN','core_WxxxQ','core_WxxxP','core_WxxxS','core_WxxxR',
'core_WxxxT','core_WxxxW','core_WxxxV','core_WxxxY','core_VxxxA','core_VxxxC','core_VxxxE','core_VxxxD','core_VxxxG','core_VxxxF','core_VxxxI','core_VxxxH','core_VxxxK','core_VxxxM',
'core_VxxxL','core_VxxxN','core_VxxxQ','core_VxxxP','core_VxxxS','core_VxxxR','core_VxxxT','core_VxxxW','core_VxxxV','core_VxxxY','core_YxxxA','core_YxxxC','core_YxxxE','core_YxxxD',
'core_YxxxG','core_YxxxF','core_YxxxI','core_YxxxH','core_YxxxK','core_YxxxM','core_YxxxL','core_YxxxN','core_YxxxQ','core_YxxxP','core_YxxxS','core_YxxxR','core_YxxxT','core_YxxxW',
'core_YxxxV','core_YxxxY','core_AxxxxA','core_AxxxxC','core_AxxxxE','core_AxxxxD','core_AxxxxG','core_AxxxxF','core_AxxxxI','core_AxxxxH','core_AxxxxK','core_AxxxxM','core_AxxxxL',
'core_AxxxxN','core_AxxxxQ','core_AxxxxP','core_AxxxxS','core_AxxxxR','core_AxxxxT','core_AxxxxW','core_AxxxxV','core_AxxxxY','core_CxxxxA','core_CxxxxC','core_CxxxxE','core_CxxxxD',
'core_CxxxxG','core_CxxxxF','core_CxxxxI','core_CxxxxH','core_CxxxxK','core_CxxxxM','core_CxxxxL','core_CxxxxN','core_CxxxxQ','core_CxxxxP','core_CxxxxS','core_CxxxxR','core_CxxxxT',
'core_CxxxxW','core_CxxxxV','core_CxxxxY','core_ExxxxA','core_ExxxxC','core_ExxxxE','core_ExxxxD','core_ExxxxG','core_ExxxxF','core_ExxxxI','core_ExxxxH','core_ExxxxK','core_ExxxxM',
'core_ExxxxL','core_ExxxxN','core_ExxxxQ','core_ExxxxP','core_ExxxxS','core_ExxxxR','core_ExxxxT','core_ExxxxW','core_ExxxxV','core_ExxxxY','core_DxxxxA','core_DxxxxC','core_DxxxxE',
'core_DxxxxD','core_DxxxxG','core_DxxxxF','core_DxxxxI','core_DxxxxH','core_DxxxxK','core_DxxxxM','core_DxxxxL','core_DxxxxN','core_DxxxxQ','core_DxxxxP','core_DxxxxS','core_DxxxxR',
'core_DxxxxT','core_DxxxxW','core_DxxxxV','core_DxxxxY','core_GxxxxA','core_GxxxxC','core_GxxxxE','core_GxxxxD','core_GxxxxG','core_GxxxxF','core_GxxxxI','core_GxxxxH','core_GxxxxK',
'core_GxxxxM','core_GxxxxL','core_GxxxxN','core_GxxxxQ','core_GxxxxP','core_GxxxxS','core_GxxxxR','core_GxxxxT','core_GxxxxW','core_GxxxxV','core_GxxxxY','core_FxxxxA','core_FxxxxC',
'core_FxxxxE','core_FxxxxD','core_FxxxxG','core_FxxxxF','core_FxxxxI','core_FxxxxH','core_FxxxxK','core_FxxxxM','core_FxxxxL','core_FxxxxN','core_FxxxxQ','core_FxxxxP','core_FxxxxS',
'core_FxxxxR','core_FxxxxT','core_FxxxxW','core_FxxxxV','core_FxxxxY','core_IxxxxA','core_IxxxxC','core_IxxxxE','core_IxxxxD','core_IxxxxG','core_IxxxxF','core_IxxxxI','core_IxxxxH',
'core_IxxxxK','core_IxxxxM','core_IxxxxL','core_IxxxxN','core_IxxxxQ','core_IxxxxP','core_IxxxxS','core_IxxxxR','core_IxxxxT','core_IxxxxW','core_IxxxxV','core_IxxxxY','core_HxxxxA',
'core_HxxxxC','core_HxxxxE','core_HxxxxD','core_HxxxxG','core_HxxxxF','core_HxxxxI','core_HxxxxH','core_HxxxxK','core_HxxxxM','core_HxxxxL','core_HxxxxN','core_HxxxxQ','core_HxxxxP',
'core_HxxxxS','core_HxxxxR','core_HxxxxT','core_HxxxxW','core_HxxxxV','core_HxxxxY','core_KxxxxA','core_KxxxxC','core_KxxxxE','core_KxxxxD','core_KxxxxG','core_KxxxxF','core_KxxxxI',
'core_KxxxxH','core_KxxxxK','core_KxxxxM','core_KxxxxL','core_KxxxxN','core_KxxxxQ','core_KxxxxP','core_KxxxxS','core_KxxxxR','core_KxxxxT','core_KxxxxW','core_KxxxxV','core_KxxxxY',
'core_MxxxxA','core_MxxxxC','core_MxxxxE','core_MxxxxD','core_MxxxxG','core_MxxxxF','core_MxxxxI','core_MxxxxH','core_MxxxxK','core_MxxxxM','core_MxxxxL','core_MxxxxN','core_MxxxxQ',
'core_MxxxxP','core_MxxxxS','core_MxxxxR','core_MxxxxT','core_MxxxxW','core_MxxxxV','core_MxxxxY','core_LxxxxA','core_LxxxxC','core_LxxxxE','core_LxxxxD','core_LxxxxG','core_LxxxxF',
'core_LxxxxI','core_LxxxxH','core_LxxxxK','core_LxxxxM','core_LxxxxL','core_LxxxxN','core_LxxxxQ','core_LxxxxP','core_LxxxxS','core_LxxxxR','core_LxxxxT','core_LxxxxW','core_LxxxxV',
'core_LxxxxY','core_NxxxxA','core_NxxxxC','core_NxxxxE','core_NxxxxD','core_NxxxxG','core_NxxxxF','core_NxxxxI','core_NxxxxH','core_NxxxxK','core_NxxxxM','core_NxxxxL','core_NxxxxN',
'core_NxxxxQ','core_NxxxxP','core_NxxxxS','core_NxxxxR','core_NxxxxT','core_NxxxxW','core_NxxxxV','core_NxxxxY','core_QxxxxA','core_QxxxxC','core_QxxxxE','core_QxxxxD','core_QxxxxG',
'core_QxxxxF','core_QxxxxI','core_QxxxxH','core_QxxxxK','core_QxxxxM','core_QxxxxL','core_QxxxxN','core_QxxxxQ','core_QxxxxP','core_QxxxxS','core_QxxxxR','core_QxxxxT','core_QxxxxW',
'core_QxxxxV','core_QxxxxY','core_PxxxxA','core_PxxxxC','core_PxxxxE','core_PxxxxD','core_PxxxxG','core_PxxxxF','core_PxxxxI','core_PxxxxH','core_PxxxxK','core_PxxxxM','core_PxxxxL',
'core_PxxxxN','core_PxxxxQ','core_PxxxxP','core_PxxxxS','core_PxxxxR','core_PxxxxT','core_PxxxxW','core_PxxxxV','core_PxxxxY','core_SxxxxA','core_SxxxxC','core_SxxxxE','core_SxxxxD',
'core_SxxxxG','core_SxxxxF','core_SxxxxI','core_SxxxxH','core_SxxxxK','core_SxxxxM','core_SxxxxL','core_SxxxxN','core_SxxxxQ','core_SxxxxP','core_SxxxxS','core_SxxxxR','core_SxxxxT',
'core_SxxxxW','core_SxxxxV','core_SxxxxY','core_RxxxxA','core_RxxxxC','core_RxxxxE','core_RxxxxD','core_RxxxxG','core_RxxxxF','core_RxxxxI','core_RxxxxH','core_RxxxxK','core_RxxxxM',
'core_RxxxxL','core_RxxxxN','core_RxxxxQ','core_RxxxxP','core_RxxxxS','core_RxxxxR','core_RxxxxT','core_RxxxxW','core_RxxxxV','core_RxxxxY','core_TxxxxA','core_TxxxxC','core_TxxxxE',
'core_TxxxxD','core_TxxxxG','core_TxxxxF','core_TxxxxI','core_TxxxxH','core_TxxxxK','core_TxxxxM','core_TxxxxL','core_TxxxxN','core_TxxxxQ','core_TxxxxP','core_TxxxxS','core_TxxxxR',
'core_TxxxxT','core_TxxxxW','core_TxxxxV','core_TxxxxY','core_WxxxxA','core_WxxxxC','core_WxxxxE','core_WxxxxD','core_WxxxxG','core_WxxxxF','core_WxxxxI','core_WxxxxH','core_WxxxxK',
'core_WxxxxM','core_WxxxxL','core_WxxxxN','core_WxxxxQ','core_WxxxxP','core_WxxxxS','core_WxxxxR','core_WxxxxT','core_WxxxxW','core_WxxxxV','core_WxxxxY','core_VxxxxA','core_VxxxxC',
'core_VxxxxE','core_VxxxxD','core_VxxxxG','core_VxxxxF','core_VxxxxI','core_VxxxxH','core_VxxxxK','core_VxxxxM','core_VxxxxL','core_VxxxxN','core_VxxxxQ','core_VxxxxP','core_VxxxxS',
'core_VxxxxR','core_VxxxxT','core_VxxxxW','core_VxxxxV','core_VxxxxY','core_YxxxxA','core_YxxxxC','core_YxxxxE','core_YxxxxD','core_YxxxxG','core_YxxxxF','core_YxxxxI','core_YxxxxH',
'core_YxxxxK','core_YxxxxM','core_YxxxxL','core_YxxxxN','core_YxxxxQ','core_YxxxxP','core_YxxxxS','core_YxxxxR','core_YxxxxT','core_YxxxxW','core_YxxxxV','core_YxxxxY','core_AxxxxxA'
,'core_AxxxxxC','core_AxxxxxE','core_AxxxxxD','core_AxxxxxG','core_AxxxxxF','core_AxxxxxI','core_AxxxxxH','core_AxxxxxK','core_AxxxxxM','core_AxxxxxL','core_AxxxxxN','core_AxxxxxQ',
'core_AxxxxxP','core_AxxxxxS','core_AxxxxxR','core_AxxxxxT','core_AxxxxxW','core_AxxxxxV','core_AxxxxxY','core_CxxxxxA','core_CxxxxxC','core_CxxxxxE','core_CxxxxxD','core_CxxxxxG',
'core_CxxxxxF','core_CxxxxxI','core_CxxxxxH','core_CxxxxxK','core_CxxxxxM','core_CxxxxxL','core_CxxxxxN','core_CxxxxxQ','core_CxxxxxP','core_CxxxxxS','core_CxxxxxR','core_CxxxxxT',
'core_CxxxxxW','core_CxxxxxV','core_CxxxxxY','core_ExxxxxA','core_ExxxxxC','core_ExxxxxE','core_ExxxxxD','core_ExxxxxG','core_ExxxxxF','core_ExxxxxI','core_ExxxxxH','core_ExxxxxK',
'core_ExxxxxM','core_ExxxxxL','core_ExxxxxN','core_ExxxxxQ','core_ExxxxxP','core_ExxxxxS','core_ExxxxxR','core_ExxxxxT','core_ExxxxxW','core_ExxxxxV','core_ExxxxxY','core_DxxxxxA',
'core_DxxxxxC','core_DxxxxxE','core_DxxxxxD','core_DxxxxxG','core_DxxxxxF','core_DxxxxxI','core_DxxxxxH','core_DxxxxxK','core_DxxxxxM','core_DxxxxxL','core_DxxxxxN','core_DxxxxxQ',
'core_DxxxxxP','core_DxxxxxS','core_DxxxxxR','core_DxxxxxT','core_DxxxxxW','core_DxxxxxV','core_DxxxxxY','core_GxxxxxA','core_GxxxxxC','core_GxxxxxE','core_GxxxxxD','core_GxxxxxG',
'core_GxxxxxF','core_GxxxxxI','core_GxxxxxH','core_GxxxxxK','core_GxxxxxM','core_GxxxxxL','core_GxxxxxN','core_GxxxxxQ','core_GxxxxxP','core_GxxxxxS','core_GxxxxxR','core_GxxxxxT',
'core_GxxxxxW','core_GxxxxxV','core_GxxxxxY','core_FxxxxxA','core_FxxxxxC','core_FxxxxxE','core_FxxxxxD','core_FxxxxxG','core_FxxxxxF','core_FxxxxxI','core_FxxxxxH','core_FxxxxxK',
'core_FxxxxxM','core_FxxxxxL','core_FxxxxxN','core_FxxxxxQ','core_FxxxxxP','core_FxxxxxS','core_FxxxxxR','core_FxxxxxT','core_FxxxxxW','core_FxxxxxV','core_FxxxxxY','core_IxxxxxA',
'core_IxxxxxC','core_IxxxxxE','core_IxxxxxD','core_IxxxxxG','core_IxxxxxF','core_IxxxxxI','core_IxxxxxH','core_IxxxxxK','core_IxxxxxM','core_IxxxxxL','core_IxxxxxN','core_IxxxxxQ',
'core_IxxxxxP','core_IxxxxxS','core_IxxxxxR','core_IxxxxxT','core_IxxxxxW','core_IxxxxxV','core_IxxxxxY','core_HxxxxxA','core_HxxxxxC','core_HxxxxxE','core_HxxxxxD','core_HxxxxxG',
'core_HxxxxxF','core_HxxxxxI','core_HxxxxxH','core_HxxxxxK','core_HxxxxxM','core_HxxxxxL','core_HxxxxxN','core_HxxxxxQ','core_HxxxxxP','core_HxxxxxS','core_HxxxxxR','core_HxxxxxT',
'core_HxxxxxW','core_HxxxxxV','core_HxxxxxY','core_KxxxxxA','core_KxxxxxC','core_KxxxxxE','core_KxxxxxD','core_KxxxxxG','core_KxxxxxF','core_KxxxxxI','core_KxxxxxH','core_KxxxxxK',
'core_KxxxxxM','core_KxxxxxL','core_KxxxxxN','core_KxxxxxQ','core_KxxxxxP','core_KxxxxxS','core_KxxxxxR','core_KxxxxxT','core_KxxxxxW','core_KxxxxxV','core_KxxxxxY','core_MxxxxxA',
'core_MxxxxxC','core_MxxxxxE','core_MxxxxxD','core_MxxxxxG','core_MxxxxxF','core_MxxxxxI','core_MxxxxxH','core_MxxxxxK','core_MxxxxxM','core_MxxxxxL','core_MxxxxxN','core_MxxxxxQ',
'core_MxxxxxP','core_MxxxxxS','core_MxxxxxR','core_MxxxxxT','core_MxxxxxW','core_MxxxxxV','core_MxxxxxY','core_LxxxxxA','core_LxxxxxC','core_LxxxxxE','core_LxxxxxD','core_LxxxxxG',
'core_LxxxxxF','core_LxxxxxI','core_LxxxxxH','core_LxxxxxK','core_LxxxxxM','core_LxxxxxL','core_LxxxxxN','core_LxxxxxQ','core_LxxxxxP','core_LxxxxxS','core_LxxxxxR','core_LxxxxxT',
'core_LxxxxxW','core_LxxxxxV','core_LxxxxxY','core_NxxxxxA','core_NxxxxxC','core_NxxxxxE','core_NxxxxxD','core_NxxxxxG','core_NxxxxxF','core_NxxxxxI','core_NxxxxxH','core_NxxxxxK',
'core_NxxxxxM','core_NxxxxxL','core_NxxxxxN','core_NxxxxxQ','core_NxxxxxP','core_NxxxxxS','core_NxxxxxR','core_NxxxxxT','core_NxxxxxW','core_NxxxxxV','core_NxxxxxY','core_QxxxxxA',
'core_QxxxxxC','core_QxxxxxE','core_QxxxxxD','core_QxxxxxG','core_QxxxxxF','core_QxxxxxI','core_QxxxxxH','core_QxxxxxK','core_QxxxxxM','core_QxxxxxL','core_QxxxxxN','core_QxxxxxQ',
'core_QxxxxxP','core_QxxxxxS','core_QxxxxxR','core_QxxxxxT','core_QxxxxxW','core_QxxxxxV','core_QxxxxxY','core_PxxxxxA','core_PxxxxxC','core_PxxxxxE','core_PxxxxxD','core_PxxxxxG',
'core_PxxxxxF','core_PxxxxxI','core_PxxxxxH','core_PxxxxxK','core_PxxxxxM','core_PxxxxxL','core_PxxxxxN','core_PxxxxxQ','core_PxxxxxP','core_PxxxxxS','core_PxxxxxR','core_PxxxxxT',
'core_PxxxxxW','core_PxxxxxV','core_PxxxxxY','core_SxxxxxA','core_SxxxxxC','core_SxxxxxE','core_SxxxxxD','core_SxxxxxG','core_SxxxxxF','core_SxxxxxI','core_SxxxxxH','core_SxxxxxK',
'core_SxxxxxM','core_SxxxxxL','core_SxxxxxN','core_SxxxxxQ','core_SxxxxxP','core_SxxxxxS','core_SxxxxxR','core_SxxxxxT','core_SxxxxxW','core_SxxxxxV','core_SxxxxxY','core_RxxxxxA',
'core_RxxxxxC','core_RxxxxxE','core_RxxxxxD','core_RxxxxxG','core_RxxxxxF','core_RxxxxxI','core_RxxxxxH','core_RxxxxxK','core_RxxxxxM','core_RxxxxxL','core_RxxxxxN','core_RxxxxxQ',
'core_RxxxxxP','core_RxxxxxS','core_RxxxxxR','core_RxxxxxT','core_RxxxxxW','core_RxxxxxV','core_RxxxxxY','core_TxxxxxA','core_TxxxxxC','core_TxxxxxE','core_TxxxxxD','core_TxxxxxG',
'core_TxxxxxF','core_TxxxxxI','core_TxxxxxH','core_TxxxxxK','core_TxxxxxM','core_TxxxxxL','core_TxxxxxN','core_TxxxxxQ','core_TxxxxxP','core_TxxxxxS','core_TxxxxxR','core_TxxxxxT',
'core_TxxxxxW','core_TxxxxxV','core_TxxxxxY','core_WxxxxxA','core_WxxxxxC','core_WxxxxxE','core_WxxxxxD','core_WxxxxxG','core_WxxxxxF','core_WxxxxxI','core_WxxxxxH','core_WxxxxxK',
'core_WxxxxxM','core_WxxxxxL','core_WxxxxxN','core_WxxxxxQ','core_WxxxxxP','core_WxxxxxS','core_WxxxxxR','core_WxxxxxT','core_WxxxxxW','core_WxxxxxV','core_WxxxxxY','core_VxxxxxA',
'core_VxxxxxC','core_VxxxxxE','core_VxxxxxD','core_VxxxxxG','core_VxxxxxF','core_VxxxxxI','core_VxxxxxH','core_VxxxxxK','core_VxxxxxM','core_VxxxxxL','core_VxxxxxN','core_VxxxxxQ',
'core_VxxxxxP','core_VxxxxxS','core_VxxxxxR','core_VxxxxxT','core_VxxxxxW','core_VxxxxxV','core_VxxxxxY','core_YxxxxxA','core_YxxxxxC','core_YxxxxxE','core_YxxxxxD','core_YxxxxxG',
'core_YxxxxxF','core_YxxxxxI','core_YxxxxxH','core_YxxxxxK','core_YxxxxxM','core_YxxxxxL','core_YxxxxxN','core_YxxxxxQ','core_YxxxxxP','core_YxxxxxS','core_YxxxxxR','core_YxxxxxT',
'core_YxxxxxW','core_YxxxxxV','core_YxxxxxY'])

fimo_file = 'fimo.tsv'
memes = leader_motif(fimo_file)
next(fh)
for line in fh:
    precursor = (line.split(',')[0] + '__'+ line.split(',')[4], line.split(',')[7].strip())
    leader, core = cut_pep(precursor, memes)
    check_pep = ''
    if "C" not in core or "X" in precursor[1] or "J" in precursor[1] or "B" in precursor[1]:
        check_pep = 'bad'
    elif "S" not in core and "T" not in core: 
        check_pep = 'bad'
    if check_pep != 'bad':
        core_params = pep_param(core)
        core_params = calc_ratios(core, core_params)
        core_params = positional_fracs(precursor[1], core_params)
        core_params = aa_pairs(core, core_params, 1)
        core_params = aa_pairs(core, core_params, 2)
        core_params = aa_pairs(core, core_params, 3)
        core_params = aa_pairs(core, core_params, 4)
        core_params = aa_pairs(core, core_params, 5)
        core_params = aa_pairs(core, core_params, 6)
        data = line.strip().split(',')
        data.extend([leader, core, len(precursor[1]), len(leader), len(core)]+ core_params)
        csv.writer(out).writerow(data)
fh.close()
out.close()

