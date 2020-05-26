import genome_overview
from genome_overview import models
import getGenomes

def test_auto_classify():
    original_classifications = {"NZ_WRXN01000057.1": "Type3", "NZ_QODI01000091.1": "Type1", "NZ_QAIK01000198.1":
        "Type2b",
                       "NZ_CP009451.1": "Type2b", "NC_015513.1": "Type3", "NC_013216.1": "Type3",
                       "NC_022663.1": "Type3", "NZ_CP004078.1": "Type3", "NZ_CP011809.2": "Type2b",
                       "NZ_VITD01000034.1": "Type3", "NZ_WSEM01000042.1": "Type3", "NZ_SOEK01000059.1": "Type2b",
                       "NZ_MQWC01000010.1": "Type3", "NZ_PVZA01000046.1": "Type2b", "NZ_BIFQ01000002.1": "Type3",
                       "NZ_QEOF01000027.1": "Type2b", "NZ_NJAK01000005.1": "Type2a", "NZ_FPBP01000034.1": "Type3",
                       "NZ_BBMZ01000055.1": "Type2b", "NZ_KI632511.1": "Type3", "NZ_CP027760.1": "Type2b",
                       "NZ_CP024793.1": "Type3", "NZ_FPJC01000063.1": "Type2b", "NZ_CP027734.1": "Type2b",
                       "NZ_SAXA01000055.1": "Type3", "NZ_NMRE01000216.1": "Type2b", "NZ_CP041692.1": "Type3",
                       "NC_008271.1": "Type3", "NZ_VCNA01000026.1": "Type3", "NZ_NBVR01000044.1": "Type1",
                       "NZ_QVIG01000004.1": "Type1", "NZ_BILZ01000167.1": "Type3", "NZ_FONE01000126.1": "Type3",
                       "NZ_WIVQ01000359.1": "Type2b", "NZ_UGTQ01000009.1": "Type1", "NZ_FOLC01000061.1": "Type1",
                       "NZ_SJOP01000106.1": "Type1", "JAAHIE010001451.1": "Type3", "NZ_KL647038.1": "Type3",
                       "NZ_WIBD01000081.1": "Type2b", "NZ_MWTQ01000158.1": "Type2b", "NZ_PHSU01000082.1": "Type2b",
                       "NZ_FUXU01000241.1": "Type1", "NZ_LT629762.1": "Type2b ", "NZ_SMSC01000216.1": "Type1",
                       "NZ_CP013429.1": "Type1", "NZ_AYLO01000184.1": "Type3", "NZ_KB905728.1": "Type3",
                       "NZNV01000041.1": "Type3", "NZ_POEF01000127.1": "Type3", "NZ_NVPT01000374.1": "Type3",
                       "NZ_KN173624.1": "Type3", "NZ_QTUF01000034.1": "Type2b", "NZ_KN266223.1": "Type3",
                       "NZ_PENW01000035.1": "Type1", "NZ_QAOO01000085.1": "Type3", "NZ_LT707062.1": "Type2b",
                       "NZ_SMKX01000348.1": "Type3", "NZ_LT855380.1": "Type1", "NZ_RQJP01000018.1": "Type3",
                       "NZ_JNYY01000082.1": "Type3", "NZ_CZQA01000015.1": "Type3", "NZ_AKJS01000210.1": "Type2b",
                       "NZ_LOYJ01000133.1": "Type2b", "NZ_BBCC01000558.1": "Type1", "NZ_OAOQ01000062.1": "Type3",
                       "MUGG01000223.1": "Type3", "NZ_CP014226.1": "Type3", "NZ_RZUM01000202.1": "Type3",
                       "NZ_CP027727.1": "Type2b", "NZ_NHML01000109.1": "Type3", "NZ_VCSG01000349.1": "Type3",
                       "NZ_AFWT01000141.1": "Type3", "NC_010830.1": "Type3", "RHGL01000187.1": "Type3",
                       "NZ_QUMQ01000001.1": "Type1", "NZ_VOBI01000057.1": "Type2b", "NZ_FNTY01000002.1": "Type2b",
                       "NZ_JYLF01000032.1": "Type2b", "NZ_CP029843.1": "Type3", "NZ_LT629732.1": "Type3",
                       "NZ_CP017687.1": "Type2b", "NZ_ONZJ01000003.1": "Type3", "NZ_AKJH01000183.1": "Type2b",
                       "NZ_SMKT01000421.1": "Type3", "NZ_RCOE01000307.1": "Type2b", "NZ_RCOD01000106.1": "Type2b",
                       "NZ_SSWO01000061.1": "Type1", "NZ_CP033700.1": "Type2b", "NZ_CP015381.1": "Type3",
                       "NZ_WIAO01000088.1": "Type3", "NZ_QARE02000057.1": "Type2b", "NZ_MASS01000135.1": "Type3",
                       "NC_020453.1": "Type1", "NC_021084.1": "Type3", "NZ_WTCR01000039.1": "Type3",
                       "NZ_CP029710.1": "Type3", "NZ_CP011521.2": "Type2b", "NZ_KI519434.1": "Type3",
                       "NZ_VHKL01000033.1": "Type3", "NZ_SNYU01000018.1": "Type2b", "LKCF01003901.1": "Type2b",
                       "NZ_JH941054.1": "Type3", "NZ_AALD02000179.1": "Type1", "NZ_KB903530.1": "Type3",
                       "NZ_FNNQ01000042.1": "Type3", "NZ_WIVR01000272.1": "Type2b", "NZ_LIUV01000061.1": "Type2b",
                       "NZ_AP020337.1": "Type2b", "NZ_KE386823.1": "Type3", "NZ_RAVY01000454.1": "Type3",
                       "NZ_RAVX01000101.1": "Type3", "NZ_CAEB01000078.1": "Type1", "NZ_MJMK01000061.1": "Type2b",
                       "NZ_BAXG01000116.1": "Type1", "NZ_RXLQ01000081.1": "Type3", "NZ_QFXK01000092.1": "Type3",
                       "MNDS01000223.1": "Type3", "NZ_LECZ01000027.1": "Type1", "NZ_CP024866.1": "Type2b",
                       "NZ_RCNX01000105.1": "Type2b", "NZ_FNBR01000014.1": "Type2b", "DLUT01000591.1": "Type1",
                       "NZ_QAIL01000624.1": "Type3", "NZ_AKJU01000280.1": "Type3", "NZ_KQ058904.1": "Type3",
                       "LADO01000066.1": "Type3", "NZ_RHHV01000044.1": "Type3", "NZ_JAABNE010000136.1": "Type1",
                       "NZ_MKMC01000100.1": "Type3", "NZ_NJAJ01000180.1": "Type2b", "NZ_FZPH01000044.1": "Type3",
                       "NZ_PPRY01000171.1": "Type2b", "NZ_CP022411.1": "Type2b", "NZ_AMZY02000039.1": "Type3",
                       "NZ_KK211074.1": "Type2b", "NZ_VEJO01000043.1": "Type2b", "NZ_CP018049.1": "Type2b",
                       "NZ_LGRC01000069.1": "Type3", "NZ_FMCR01000011.1": "Type3", "DPJM01000861.1": "Type3",
                       "NC_008027.1": "Type2b", "NC_017448.1": "Type3", "NZ_AZNV01000101.1": "Type3",
                       "NZ_CP028042.1": "Type2b", "NZ_AP017313.1": "Type3", "NZ_QUMR01000030.1": "Type2b",
                       "NZ_FODL01000038.1": "Type2b", "QJPH01000597.1": "Type3", "NZ_FORG01000094.1": "Type2b",
                       "NZ_CP007231.1": "Type2b", "NZ_NEVM01000005.1": "Type1", "NZ_CP009289.1": "Type3",
                       "NZ_VDFY01000321.1": "Type3", "NZ_BJMN01000147.1": "Type3", "NZ_NIBS01000173.1": "Type2a",
                       "JMHR01000441.1": "Type3", "JAABOU010001898.1": "Type3", "NZ_VIUK01000070.1": "Type1",
                       "JAAHFV010000687.1": "Type3", "NZ_PVZG01000092.1": "Type3", "NZ_KI911557.1": "Type3",
                       "NZ_JOIX01000228.1": "Type3", "JEMY01000085.1": "Type3", "NZ_KI421497.1": "Type3",
                       "NZ_CP045011.1": "Type1", "NZ_SSNI01000125.1": "Type3", "NZ_SZWE01000003.1": "Type3",
                       "NZ_CDSC02000515.1": "Type3", "NZ_FMYF01000036.1": "Type3", "NZ_RCWL01000032.1": "Type3",
                       "NZ_PHHE01000001.1": "Type2b", "NZ_LKBY01000178.1": "Type3", "NZ_AP018150.1": "Type2b",
                       "NZ_FMXV01000075.1": "Type2b", "NZ_KB897775.1": "Type3", "NZ_PENZ01000052.1": "Type1",
                       "NZ_NIRH01000051.1": "Type1", "NZ_PENX01000027.1": "Type3", "NZ_CP047651.1": "Type2b",
                       "NZ_SNXZ01000022.1": "Type3", "NC_017565.1": "Type1", "NZ_AP018449.1": "Type3",
                       "NZ_MSSW01000136.1": "Type3", "NZ_LR134373.1": "Type2b", "NZ_CP038274.1": "Type3",
                       "NZ_ASSC01000896.1": "Type3", "NZ_FOSU01000047.1": "Type3", "NZ_LAIJ01000019.1": "Type3",
                       "NZ_FUYT01000034.1": "Type2b", "NZ_MBLO01000280.1": "Type3", "NZ_KE384514.1": "Type3",
                       "NZ_CP009747.1": "Type2b", "NZ_QGSY01000345.1": "Type3", "NZ_QGSZ01000408.1": "Type3",
                       "NZ_JH725405.1": "Type3", "NZ_OUNR01000022.1": "Type3", "NZ_CP005927.1": "Type1",
                       "NZ_CP043925.1": "Type1", "NZ_ASRX01000182.1": "Type1", "NZ_BAHC01000261.1": "Type3",
                       "NZ_PVTJ01000022.1": "Type3", "NZ_LR590468.1": "Type3", "NZ_LLWH01000241.1": "Type2b",
                       "NC_017447.1": "Type1", "NZ_MWPQ01000095.1": "Type3", "NZ_RHQN01000027.1": "Type2b",
                       "NZ_QAOQ01000022.1": "Type3", "NZ_VUAZ01000259.1": "Type2b", "NZ_CP033931.1": "Type3",
                       "NZ_LS999839.1": "Type3", "NZ_KI421431.1": "Type3", "NZ_CP027759.1": "Type2b",
                       "NZ_QTUH01000043.1": "Type2b", "NZ_BCBA01000109.1": "Type2b", "NC_015559.1": "Type3",
                       "NZ_AKJQ01000091.1": "Type2b", "NZ_QTPO01000204.1": "Type3", "NZ_PDUD01000143.1": "Type3",
                       "JAAAKW010000055.1": "Type2b", "NZ_WSTC01000100.1": "Type3", "NZ_SOCG01000010.1": "Type2b",
                       "NZ_QTPW01000119.1": "Type3", "NZ_JAABMA010000050.1": "Type1", "NZ_LT629778.1": "Type2b",
                       "NZ_PJBP01000186.1": "Type3", "NZ_QFRW01000331.1": "Type3", "NZ_MKQR01000032.1": "Type3",
                       "NZ_CP046054.1": "Type3", "NZ_UPHT01000230.1": "Type3", "NZ_MCHY01000013.1": "Type3",
                       "NZ_MUNY01000094.1": "Type3", "NZ_RAWG01000802.1": "Type3", "NZ_JYLO01000042.1": "Type2b",
                       "NZ_WSQA01000032.1": "Type3", "NZ_CP028923.1": "Type3", "NZ_MUBJ01000149.1": "Type2a",
                       "NZ_JXRA01000201.1": "Type2b", "NZ_JYLD01000037.1": "Type3", "NZ_BDBY01000492.1": "Type3",
                       "NZ_LIPP01000561.1": "Type3", "NZ_AHAM01000375.1": "Type3", "NZ_BILY01000094.1": "Type3",
                       "NZ_SODV01000004.1": "Type3", "NZ_FYEA01000033.1": "Type2b", "NZ_WMBA01000144.1": "Type3",
                       "NZ_FNCO01000054.1": "Type2b", "NZ_FMUL01000047.1": "Type2b", "NZ_FCON02000657.1": "Type3",
                       "NZ_CP023969.1": "Type2b", "NZ_JJML01000118.1": "Type3", "NZ_JAABLV010000025.1": "Type1",
                       "NZ_JAABLU010000039.1": "Type1", "NZ_FORB01000034.1": "Type3", "NZ_JYLH01000043.1": "Type2b",
                       "NZ_PGGO01000087.1": "Type3", "NZ_LMGQ01000029.1": "Type2b", "NZ_JAABNK010000025.1": "Type1",
                       "NZ_KZ679081.1": "Type3", "NKIG01000124.1": "Type3", "NZ_LMGK01000026.1": "Type2b",
                       "WASQ01000153.1": "Type3", "NZ_BAOS01000047.1": "Type3", "NZ_BCQP01000133.1": "Type3",
                       "NZ_CP010898.2": "Type2b", "NC_021184.1": "Type3", "NZ_FOZR01000085.1": "Type3",
                       "NC_009253.1": "Type3", "NZ_QKLY01000024.1": "Type2b", "NZ_LVYD01000134.1": "Type3",
                       "NZ_VFIO01000040.1": "Type2b", "QQTZ01000066.1": "Type3", "NZ_FNJL01000093.1": "Type3",
                       "NZ_MVHE01000525.1": "Type3", "NZ_FMWY01000062.1": "Type3", "NZ_CP010408.1": "Type3",
                       "NZ_LT605205.1": "Type3", "LKBL01002861.1": "Type1", "NZ_KB944506.1": "Type3",
                       "NZ_CP029618.1": "Type3", "NZ_FPBO01000103.1": "Type3", "NZ_QJUG01000493.1": "Type3",
                       "NZ_QAJM01000070.1": "Type2b", "LGGF01000107.1": "Type3", "NZ_WUNA01000095.1": "Type3",
                       "NZ_MKCS01000005.1": "Type3", "NZ_CM002331.1": "Type2b", "NZ_CP023695.1": "Type3",
                       "NZ_SMFY01000011.1": "Type3", "NZ_QSNX01000060.1": "Type2b", "NZ_LMXH01000018.1": "Type3",
                       "NZ_CP014135.1": "Type2b", "NZ_JXDG01000126.1": "Type2b", "NZ_PIQI01000031.1": "Type1",
                       "NZ_QKTW01000033.1": "Type3", "NZ_LOWA01000060.1": "Type2b", "NZ_CP042382.1": "Type3",
                       "NZ_CP013461.1": "Type3", "NZ_LWBP01000264.1": "Type3", "NZ_LYRP01000050.1": "Type1",
                       "NZ_SEIT01000119.1": "Type2b", "NZ_JYLN01000037.1": "Type2b", "NZ_QTTH01000050.1": "Type2b",
                       "NZ_NITZ01000138.1": "Type1", "NZ_VOLC01000068.1": "Type3", "NZ_LFCV01000266.1": "Type1",
                       "NZ_MVIF01000387.1": "Type3", "NZ_VRLS01000101.1": "Type3", "NZ_MULM01000185.1": "Type3",
                       "JAABLX010000056.1": "Type1", "NC_016905.1": "Type1", "NZ_SOZA01000107.1": "Type2b",
                       "NZ_MJML01000057.1": "Type2b", "NZ_WTYM01000063.1": "Type3", "NZ_QOIO01000184.1": "Type3",
                       "NZ_FQUQ01000022.1": "Type3", "NZ_FAOZ01000083.1": "Type3", "NZ_JNZS01000088.1": "Type3",
                       "NZ_KQ257877.1": "Type3", "NZ_KB892704.1": "Type3", "AP018273.1": "Type3",
                       "NZ_LT985385.1": "Type3", "NZ_PYAC01000043.1": "Type3", "NZ_FOBB01000023.1": "Type3",
                       "NZ_JPMW01000007.1": "Type3", "NZ_CPYD01000045.1": "Type2a", "NZ_CP021135.1": "Type2b",
                       "MNDA01000671.1": "Type3", "NZ_QXQA01000045.1": "Type3", "NZ_OCSV01000008.1": "Type3",
                       "NZ_FXWP01000029.1": "Type1", "NZ_AWXZ01000044.1": "Type3", "NZ_UYJA01000022.1": "Type2b",
                       "NZ_LNTU01000041.1": "Type3", "NZ_QNVV01000061.1": "Type3", "NZ_CP022121.1": "Type3",
                       "NZ_SAIQ01000015.1": "Type3", "NZ_VCRA01000094.1": "Type3", "NZ_CP029197.1": "Type3",
                       "NZ_PODL01000171.1": "Type2b", "NZ_FOAF01000023.1": "Type3", "NZ_QKWJ01000248.1": "Type3",
                       "NZ_CP029608.1": "Type2b", "NZ_JFHN01000075.1": "Type1", "NZ_FXBM01000005.1": "Type3",
                       "NZ_CP048209.1": "Type3", "NZ_VJZE01000939.1": "Type3", "NC_013131.1": "Type3",
                       "NZ_JH651384.1": "Type3", "NZ_PYAW01000034.1": "Type3", "NZ_WMJZ01000174.1": "Type1",
                       "NZ_SNZP01000034.1": "Type3", "NZ_CP010896.1": "Type2b", "NZ_SMJU01000044.1": "Type3",
                       "NZ_FAOS01000004.1": "Type3", "NZ_RHLK01000044.1": "Type3", "NZ_VSFF01000027.1": "Type3",
                       "NZ_RQPI01000039.1": "Type3", "NZ_FSRS01000002.1": "Type3", "NZ_QLTF01000036.1": "Type2b",
                       "NZ_FNKR01000003.1": "Type3", "NZ_SMSL01000028.1": "Type3", "NZ_VZZK01000111.1": "Type3",
                       "NZ_CP028272.1": "Type1", "NZ_VDCQ01000167.1": "Type3", "NZ_OGTP01000072.1": "Type3",
                       "NZ_MPIN01000042.1": "Type3", "NZ_CDPK01000072.1": "Type1", "NZ_CP026364.1": "Type1",
                       "NZ_LXEN01000293.1": "Type3", "NZ_CABPSP010000048.1": "Type2b", "NZ_CP019686.1": "Type1",
                       "NZ_SJSL01000015.1": "Type3", "CABPSQ010000036.1": "Type2b", "JAAHFO010001461.1": "Type3",
                       "NZ_SOCQ01000042.1": "Type2b", "MNJJ01000259.1": "Type3", "NZ_CP012159.1": "Type3",
                       "NZ_QLTJ01000050.1": "Type2b", "NZ_JNWO01000211.1": "Type3", "NZ_CP013341.1": "Type3",
                       "NC_017807.1": "Type2b", "NZ_PYBV01000203.1": "Type3", "NZ_KE332397.1": "Type3",
                       "NZ_RCSU01000043.1": "Type3", "NZ_QLLL01000025.1": "Type3", "NZ_QTUB01000001.1": "Type2a",
                       "NZ_JAAGLX010000268.1": "Type3", "NZ_FODH01000041.1": "Type3", "NZ_FNVU01000044.1": "Type3",
                       "NZ_FXAS01000142.1": "Type1", "NZ_CP038255.1": "Type3", "NZ_QVNU01000027.1": "Type3",
                       "NZ_VOIW01000021.1": "Type2b", "NZ_LT629795.1": "Type1", "NZ_CP026110.1": "Type3",
                       "NZ_AEDD01000040.1": "Type3", "NZ_FNAD01000036.1": "Type3", "NZ_PVZV01000026.1": "Type3",
                       "NZ_PVNL01000192.1": "Type3", "NZ_LXYR01000213.1": "Type3", "NZ_LN623556.1": "Type2b",
                       "NZ_FNGF01000016.1": "Type3", "NZ_CP022961.1": "Type3", "NZ_CXPG01000027.1": "Type3",
                       "NZ_CP017236.1": "Type1", "NZ_CP012332.1": "Type3", "NZ_NIRS01000013.1": "Type2b",
                       "NZ_FMDM01000037.1": "Type3", "NZ_PTJB01000052.1": "Type3", "NZ_JAAFZB010000096.1": "Type3",
                       "NZ_CP016211.1": "Type3", "NZ_PQKR01000051.1": "Type2b", "NZ_NCXP01000142.1": "Type3",
                       "NZ_ANMG01000154.1": "Type3", "NC_020209.1": "Type2b", "NZ_JZSQ01000140.1": "Type3",
                       "NZ_LT629705.1": "Type2b", "NZ_PYAL01000013.1": "Type2b", "MKSF01000039.1": "Type3",
                       "LAQJ01000315.1": "Type3", "NZ_SMKU01000805.1": "Type3", "NZ_LT828648.1": "Type3",
                       "NZ_CP007215.2": "Type3", "NZ_WNKZ01000359.1": "Type3", "NZ_LR590482.1": "Type2b",
                       "NZ_LT907981.1": "Type3", "NZ_QAIP01000588.1": "Type1", "NZ_LNCD01000152.1": "Type3",
                       "NZ_KE384562.1": "Type3", "NZ_ATXB01000005.1": "Type3", "NZ_SMKK01000563.1": "Type3",
                       "NC_019762.1": "Type3", "NZ_JOGP01000180.1": "Type3", "KZ266893.1": "Type3",
                       "NZ_FNON01000025.1": "Type3", "NZ_SHKK01000001.1": "Type3", "NZ_FNUD01000002.1": "Type1",
                       "NZ_FQYP01000028.1": "Type3", "NZ_QGTQ01000078.1": "Type1", "NZ_JFJW01000247.1": "Type1",
                       "NZ_FOVS01000095.1": "Type2b", "NZ_CP012540.1": "Type3", "NZ_JUHO01000001.1": "Type2b",
                       "DNUG01000139.1": "Type3", "NZ_CP038630.1": "Type3", "NC_013954.1": "Type1",
                       "NZ_FXYF01000056.1": "Type3"}



    reclassify = False


    if reclassify:
        print ("Deleting all tags")
        delete_all_tags()
        all_genomes = models.GenomeRecords.objects()
        getGenomes.tag_as_simple(all_genomes, "hidden")
        queries = models.GenomeRecords.objects(tags="Simple").timeout(False)
        print ("Classifying the genomes")
        genome_overview.classify_genomes(queries)

    else:

        queries = models.GenomeRecords.objects(tags="Simple").timeout(False)


    count = 0
    diff_count = 0
    for query in queries:
        if query.name in original_classifications:
            count += 1

            new_classification = models.GenomeTags.objects().get(tag_id=query.name)


            if original_classifications[query.name] != new_classification.tags[0].split("Auto_")[1]:
                diff_count += 1
                print ("DIFFERENT")
                print(query.name)
                print("Original was " + original_classifications[query.name])
                print("Automatic classification was " + new_classification.tags[0].split("Auto_")[1])

    print ("Wrong: " + str(diff_count))
    print ("Correct " + str(count - diff_count))
    print ("Total " + str(count))




def delete_all_tags():
    queries = models.GenomeRecords.objects().all()

    queries.update(tags=[])

    for query in queries:
        print (query.name)

        for hit in query.hits:
            hit.tags = []

        query.save()

    models.GenomeTags.objects().all().delete()