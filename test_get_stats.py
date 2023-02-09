#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:17:50 2023

@author: prowe
"""

import numpy as np

from get_inversion import get_stats


class TestGetStats:
    def test_one_inv(self):

        # Setup
        # fmt: off
        height = np.array(
            [0.01510407, 0.03559937, 0.05791085, 0.08242772, 0.10915002,
               0.13846698, 0.17037867, 0.20540406, 0.24367298, 0.28518555,
               0.33072028, 0.38027733, 0.4342461 , 0.49301601, 0.55697655,
               0.62651726, 0.70189799, 0.78376789, 0.87251673, 0.96840458,
               1.07182137, 1.18341666, 1.30345084, 1.43231419, 1.57039713,
               1.7180902 , 1.87552454, 2.04348039, 2.2216998 , 2.41044411,
               2.60971521, 2.81899574, 3.03854735, 3.26772293, 3.50639468,
               3.75365565, 4.00911821, 4.27226477, 4.54231781, 4.81862954,
               5.10029211, 5.38561778, 5.67252836, 5.95985465, 6.24746679,
               6.53627474, 6.82653879, 7.11760932])

        temp = np. array([
            243.67555647, 245.63861233, 246.53219989, 247.12225136,
               248.71463385, 250.68620007, 252.68897092, 254.53666415,
               255.91155972, 256.96968087, 258.43629919, 260.09203657,
               261.56149167, 263.50185324, 265.08288859, 266.73389799,
               268.19957072, 269.56028236, 270.87465983, 271.99802706,
               272.86513795, 273.25283203, 273.83342754, 274.77902285,
               275.63762339, 276.24374998, 276.8035424 , 277.32551101,
               277.75859366, 278.14439654, 278.43942228, 278.63326931,
               278.74484956, 278.76659825, 278.79307492, 278.86399457,
               278.96895565, 279.10039339, 279.24790626, 279.38879996,
               279.52874807, 279.66018582, 279.77554844, 279.87578155,
               279.96088512, 280.01289286, 280.01951203, 279.93535405])
        # fmt: on

        # Run
        invdepth, invstrength = get_stats(temp, height, 1)

        # Assert
        assert np.isclose(invdepth, 6.81143472)
        assert np.isclose(invstrength, 36.34395555999998)

    def test_mini_noninv_near_surf(self):

        # Setup
        # fmt: off
        temp = np.array(
            [233.23098273, 233.25193702, 233.17631935, 233.08885794,
             232.96859851, 232.96131006, 233.44143674, 234.7606463,
             236.78319134, 239.01163511, 241.13804057, 243.17242935,
             245.15579897, 246.98611113, 248.81095694, 250.57476199,
             252.23197345, 253.74979328, 255.15190897, 256.4820512,
             257.70833301, 258.84442025, 259.89031291, 260.86605424,
             261.83632922, 262.74556344, 263.56551413, 264.33900095,
             265.04780277, 265.70194121, 266.27499564, 266.79338669,
             267.28353499, 267.73906315, 268.14357216, 268.49706202,
             268.82777546, 269.16942158, 269.52928883, 269.88824502,
             270.22442481, 270.54329452, 270.84485417, 271.12090423,
             271.376     , 271.61925204, 271.86705936, 272.11759985,
             271.03]
            )
        height = np.array(
            [0.01341739, 0.03326936, 0.05494994, 0.07858979, 0.10445015,
             0.13266169, 0.16335509, 0.19692225, 0.23349388, 0.27333131,
             0.31682651, 0.36411026, 0.41557459, 0.47148095, 0.53235206,
             0.59831886, 0.66977355, 0.74723903, 0.83110767, 0.9215106,
             1.01897098, 1.12388139, 1.23663454, 1.35762324, 1.4869791,
             1.62522581, 1.77236454, 1.92904993, 2.09489133, 2.27028225,
             2.45522426, 2.64945762, 2.85272257, 3.06489011, 3.28543906,
             3.51384818, 3.74959612, 3.99189989, 4.23984555, 4.49251888,
             4.74913623, 5.00878287, 5.27080543, 5.53507368, 5.80158816,
             6.07048029, 6.3416198 , 6.61513815, 6.8]
            )
        # fmt: on

        # Run
        invdepth, invstrength = get_stats(temp, height)

        # Assert
        assert np.isclose(invdepth, 6.60172076)
        assert np.isclose(invstrength, 38.88661711999998)

    def test_mini_noninv_near_surf_inv_to_top(self):

        # fmt: off
        temp = np.array(
            [233.23098273, 233.25193702, 233.17631935, 233.08885794,
             232.96859851, 232.96131006, 233.44143674, 234.7606463,
             236.78319134, 239.01163511, 241.13804057, 243.17242935,
             245.15579897, 246.98611113, 248.81095694, 250.57476199,
             252.23197345, 253.74979328, 255.15190897, 256.4820512,
             257.70833301, 258.84442025, 259.89031291, 260.86605424,
             261.83632922, 262.74556344, 263.56551413, 264.33900095,
             265.04780277, 265.70194121, 266.27499564, 266.79338669,
             267.28353499, 267.73906315, 268.14357216, 268.49706202,
             268.82777546, 269.16942158, 269.52928883, 269.88824502,
             270.22442481, 270.54329452, 270.84485417, 271.12090423,
             271.376     , 271.61925204, 271.86705936, 272.11759985]
            )
        height = np.array(
            [0.01341739, 0.03326936, 0.05494994, 0.07858979, 0.10445015,
             0.13266169, 0.16335509, 0.19692225, 0.23349388, 0.27333131,
             0.31682651, 0.36411026, 0.41557459, 0.47148095, 0.53235206,
             0.59831886, 0.66977355, 0.74723903, 0.83110767, 0.9215106,
             1.01897098, 1.12388139, 1.23663454, 1.35762324, 1.4869791,
             1.62522581, 1.77236454, 1.92904993, 2.09489133, 2.27028225,
             2.45522426, 2.64945762, 2.85272257, 3.06489011, 3.28543906,
             3.51384818, 3.74959612, 3.99189989, 4.23984555, 4.49251888,
             4.74913623, 5.00878287, 5.27080543, 5.53507368, 5.80158816,
             6.07048029, 6.3416198 , 6.61513815]
            )
        # fmt: on

        # Run
        invdepth, invstrength = get_stats(temp, height)

        # Assert
        assert np.isclose(invdepth, 6.60172076)
        assert np.isclose(invstrength, 38.88661711999998)

    def test_several_switches(self):
        # fmt: off
        height = np.array(
            [0.01082343, 0.03131869, 0.05363015, 0.07814698, 0.10473953,
               0.13405645, 0.16596809, 0.20099343, 0.2392623 , 0.28090455,
               0.32643921, 0.37599619, 0.43022435, 0.48912393, 0.55334387,
               0.62301424, 0.69865437, 0.78078369, 0.86953244, 0.96542021,
               1.0688369 , 1.18056185, 1.30072569, 1.42971871, 1.56806109,
               1.71588384, 1.87370745, 2.041793  , 2.22027193, 2.40901616,
               2.60802752, 2.81730794, 3.03646992, 3.26564534, 3.50431693,
               3.75196735, 4.00820905, 4.27200499, 4.54309722, 4.82044832,
               5.10328036, 5.38951584, 5.67733641, 5.96544289, 6.25383536,
               6.54290385, 6.8334285 , 7.12475967])
    
        temp = np.array([
            243.67744766, 245.59700613, 246.24568451, 247.27260101,
               249.05410257, 251.08429569, 253.12678155, 255.09645657,
               256.63115775, 257.67793175, 258.89491291, 260.41070219,
               261.88299408, 263.48294133, 264.93632132, 266.22138534,
               267.56413067, 269.04587851, 270.56355898, 271.87793645,
               272.88972343, 273.64336289, 274.38281842, 275.17995526,
               275.90806364, 276.52553738, 277.05790754, 277.50139173,
               277.91650807, 278.33918917, 278.78267337, 279.18266019,
               279.51550973, 279.72164951, 279.70651999, 279.53063926,
               279.36610568, 279.28289329, 279.27627412, 279.30464198,
               279.34530258, 279.39825592, 279.44364449, 279.47957711,
               279.50038021, 279.50416259, 279.46444759, 279.30747877])
        # fmt: on

        # Run
        invdepth, invstrength = get_stats(temp, height)

        # Assert
        assert np.isclose(invdepth, 3.25482191)
        assert np.isclose(invstrength, 36.04420185000001)

    def test_inv_to_top_with_thick_noninv(self):
        # fmt: off
        height = np.array(
            [0.00861825, 0.02898378, 0.0511655 , 0.0755526 , 0.10214512,
               0.13120257, 0.16298447, 0.19775033, 0.23575971, 0.27727218,
               0.32267706, 0.37223399, 0.42646208, 0.48549132, 0.54984093,
               0.61964096, 0.69541076, 0.77754   , 0.86667791, 0.9628251 ,
               1.06676075, 1.17874515, 1.29929823, 1.42868051, 1.56728241,
               1.71562427, 1.87383724, 2.0423122 , 2.22105076, 2.41044411,
               2.61062397, 2.82146253, 3.04257239, 3.27382584, 3.51470566,
               3.76495427, 4.02418426, 4.29174835, 4.56660942, 4.84760014,
               5.13342288, 5.42277956, 5.7137219 , 6.00482077, 6.29568628,
               6.58631838, 6.876457  , 7.166102  ])

        temp = np.array(
            [241.37114071, 243.75876886, 246.16247212, 248.55482825,
               250.8809927 , 253.11354321, 255.20425443, 257.104901  ,
               258.85046993, 260.56199743, 262.06927635, 263.30516941,
               264.42191747, 265.58121731, 266.8738461 , 268.23361215,
               269.55082641, 270.67041125, 271.56210762, 272.32520303,
               273.11004714, 273.94784458, 274.85088809, 275.6953047 ,
               276.36384058, 276.91606824, 277.43236328, 277.88530343,
               278.27867107, 278.62286777, 278.91411112, 279.12592447,
               279.35381294, 279.58548379, 279.7632557 , 279.64127391,
               278.96895565, 278.69000503, 278.20586024, 277.61297198,
               277.4881534 , 277.5477259 , 277.58460412, 277.74157294,
               277.8314045 , 277.94014796, 278.01768677, 278.02336034])
        # fmt: on

        # Run
        invdepth, invstrength = get_stats(temp, height)

        # Assert
        assert np.isclose(invdepth, 3.50608741)
        assert np.isclose(invstrength, 38.39211499)

    def test_simple_inv(self):
        # fmt: off
        height = np.array(
            [6.41307063e-03, 2.65191571e-02, 4.83117038e-02, 7.24393403e-02,
               9.85129594e-02, 1.27440657e-01, 1.58833345e-01, 1.93209989e-01,
               2.30700414e-01, 2.71693909e-01, 3.16450068e-01, 3.65358242e-01,
               4.18548357e-01, 4.76669302e-01, 5.39980825e-01, 6.08872458e-01,
               6.83733793e-01, 7.65084242e-01, 8.53183803e-01, 9.48552054e-01,
               1.05144917e+00, 1.16239492e+00, 1.28164970e+00, 1.40986330e+00,
               1.54742615e+00, 1.69459902e+00, 1.85190240e+00, 2.01985701e+00,
               2.19833470e+00, 2.38759689e+00, 2.58712620e+00, 2.79718423e+00,
               3.01751337e+00, 3.24811576e+00, 3.48873391e+00, 3.73872075e+00,
               3.99742908e+00, 4.26447137e+00, 4.53894039e+00, 4.81992867e+00,
               5.10639854e+00, 5.39666232e+00, 5.68890175e+00, 5.98220767e+00,
               6.27606047e+00, 6.57059029e+00, 6.86540733e+00, 7.16038164e+00])

        temp = np. array(
            [245.73789983, 247.7387795 , 249.52879141, 251.2024951 ,
               252.71733878, 254.17828353, 255.78863233, 257.38952518,
               258.77293112, 259.92277501, 261.06410854, 262.2612322 ,
               263.48672371, 264.66588106, 265.77411876, 266.80292645,
               267.86388438, 268.99387077, 270.14087788, 271.15455604,
               271.87415407, 272.13986635, 272.16161504, 272.6334671 ,
               273.37197704, 273.61499503, 273.87503374, 274.53789605,
               275.16388014, 275.44093956, 275.67071922, 275.59696279,
               275.38325825, 275.11470918, 274.83386738, 274.67973534,
               274.4877795 , 274.39132877, 274.3667433 , 274.38281842,
               274.38754639, 274.39321996, 274.39321996, 274.3856552 ,
               274.36768889, 274.36485211, 274.40078473, 274.41496866])
        # fmt: on

        # Run
        invdepth, invstrength = get_stats(temp, height)

        # Assert
        assert np.isclose(invdepth, 2.5807131293700003)
        assert np.isclose(invstrength, 29.93281939000002)

    def test_all_inv(self):
        # fmt: off
        height = np.array(
            [0.00939655, 0.02937293, 0.05142494, 0.07542288, 0.10175596,
               0.13042425, 0.16194669, 0.19619364, 0.23355437, 0.27441818,
               0.31891492, 0.36743394, 0.42010515, 0.47744771, 0.53972135,
               0.60731557, 0.68061994, 0.75989436, 0.84578802, 0.93856096,
               1.03860307, 1.14643407, 1.26244406, 1.38689348, 1.52043222,
               1.66319118, 1.81595037, 1.9787112 , 2.15212425, 2.33606142,
               2.53065437, 2.73577522, 2.95142601, 3.17747902, 3.41367676,
               3.65885269, 3.91196989, 4.17225093, 4.43865838, 4.71028443,
               4.98609105, 5.26478002, 5.54479278, 5.82534997, 6.10619179,
               6.38718834, 6.66833966, 6.94964577])
    
        temp = np.array(
            [234.51368356, 236.44080679, 238.34239895, 240.16266991,
               241.92336837, 243.74080255, 245.66225221, 247.72837795,
               249.94769013, 252.21806446, 254.46290771, 256.65006966,
               258.73037933, 260.40597421, 261.37237261, 261.87164693,
               262.48722948, 263.26450882, 263.84510433, 264.10230626,
               264.27724139, 264.4758164 , 264.7982644 , 265.3835879 ,
               265.42519409, 266.3906469 , 266.95043932, 267.42229137,
               267.95655272, 268.49837883, 269.09221268, 269.71630558,
               270.25434931, 270.83116245, 271.35407665, 271.75406347,
               272.18809171, 272.65048782, 273.07316892, 273.46086299,
               273.80978766, 274.15398435, 274.46792199, 274.75349178,
               275.00785691, 275.24898372, 275.46930742, 275.68112077])
        # fmt: on

        # Run
        invdepth, invstrength = get_stats(temp, height)

        # Assert
        assert np.isclose(invdepth, 6.94024922)
        assert np.isclose(invstrength, 41.16743721)

    def test_inv_goes_to_top(self):
        # fmt: off
        height = np.array([
            3.31522919, 3.33315759, 3.35251518, 3.37304214, 3.39460859,
               3.4176043 , 3.44215923, 3.46840333, 3.49646657, 3.52673876,
               3.5593499 , 3.59442999, 3.6319791 , 3.67264697, 3.71643374,
               3.76359938, 3.81440392, 3.86910742, 3.92796995, 3.99125163,
               4.05934256, 4.13224306, 4.2103433 , 4.29403354, 4.38357413,
               4.47909554, 4.58072826, 4.68886279, 4.80349981, 4.92464   ,
               5.05254408, 5.18708289, 5.32825724, 5.47580803, 5.62960612,
               5.78939239, 5.95503772, 6.12602289, 6.30195865, 6.48258573,
               6.66712462, 6.85492569, 7.04481891, 7.23628436, 7.42893205,
               7.6228923 , 7.81842555, 8.01540197])
    
        temp = np.array([223.36191688, 224.45054535, 225.60226907, 226.81708803,
               228.031907  , 229.17986383, 230.22988177, 231.20079523,
               232.11897235, 233.01454819, 233.93743392, 234.88857125,
               235.83876685, 236.80873859, 237.74575011, 238.63944251,
               239.47945687, 240.26202629, 240.93064758, 241.50980546,
               241.99290789, 242.38843034, 242.68036358, 242.90825985,
               243.11732172, 243.34521799, 243.54957126, 243.78782645,
               244.038324  , 244.30859768, 244.53366879, 244.66268445,
               244.68528573, 244.68716917, 244.74555582, 244.92542436,
               245.19193116, 245.42076915, 245.72965335, 246.03100379,
               246.30504435, 246.56213395, 246.8079229 , 247.03770261,
               247.25147308, 247.45017603, 247.64699553, 247.83628128])
        # fmt: on

        # Run
        invdepth, invstrength = get_stats(temp, height)

        # Assert
        assert np.isclose(invdepth, 4.700172779999999)
        assert np.isclose(invstrength, 24.474364400000013)

    def test_one_non_inv_level(self):
        # fmt: off
        height = np.array([
            0.00988003, 0.0298851 , 0.0518263 , 0.07570365, 0.10177534,
               0.13029957, 0.16127639, 0.1950931 , 0.23200794, 0.27227916,
               0.31603596, 0.36392389, 0.41581405, 0.47222296, 0.53353814,
               0.60001803, 0.67205023, 0.75015151, 0.8345805 , 0.92572501,
               1.02397294, 1.12984137, 1.24371843, 1.36586319, 1.496664  ,
               1.63650931, 1.78565857, 1.9443714 , 2.11329495, 2.29243087,
               2.48229752, 2.68276766, 2.89358492, 3.11436379, 3.34484793,
               3.58465183, 3.83313148, 4.08964277, 4.35354142, 4.62327826,
               4.89743287, 5.17458437, 5.45305283, 5.7327093 , 6.01394201,
               6.29765662, 6.58359512, 6.87214613])
        
        temp = np.array([
            245.68374416, 245.51908829, 245.55281299, 245.44469558,
               245.44568749, 246.40485748, 248.25971568, 250.35957388,
               252.41380454, 254.20716477, 255.4281971 , 256.34074766,
               257.15311604, 258.04483664, 258.97722526, 259.74495801,
               260.34109159, 260.68627376, 260.83307537, 260.97789317,
               261.31514012, 261.92317653, 262.65420888, 263.4100388 ,
               264.1926501 , 264.93062577, 265.6199982 , 266.19827752,
               266.81424915, 267.48378353, 268.10669848, 268.70679966,
               269.11943122, 269.66101014, 270.09347976, 270.39204249,
               270.81756879, 271.27979548, 271.70333797, 272.09712338,
               272.45916789, 272.79839865, 273.11382373, 273.39155651,
               273.65242694, 273.89841883, 274.13052408, 274.34080747])
        # fmt: on

        # Run
        invdepth, invstrength = get_stats(temp, height)

        # Assert
        assert np.isclose(invdepth, 6.8622661)
        assert np.isclose(invstrength, 28.657063310000012)

    def test_invstrength_lt0(self):
        # fmt: off
        height = np.array([
            3.52472066, 3.54105616, 3.55855858, 3.57722793, 3.59706423,
               3.61806751, 3.64023779, 3.66474196, 3.6904132 , 3.71841843,
               3.7487577 , 3.78143108, 3.81643865, 3.85494742, 3.89695752,
               3.94130211, 3.98914829, 4.0404962 , 4.09651305, 4.15603204,
               4.21905335, 4.28791139, 4.36143942, 4.44080491, 4.52484111,
               4.6135484 , 4.70926166, 4.80964696, 4.91703937, 5.02910489,
               5.14817882, 5.27309448, 5.40385258, 5.54162137, 5.68289895,
               5.83118877, 5.98415636, 6.1418024 , 6.30295975, 6.46879683,
               6.63814642, 6.80984105, 6.9838811 , 7.1590988 , 7.33549436,
               7.51306795, 7.69181979, 7.87291847])

        temp = np.array([
            229.3011468 , 242.09456617, 247.47562985, 248.28665853,
               247.14377775, 244.99790003, 243.91603242, 244.83718242,
               240.49631883, 234.63459045, 231.70744657, 231.73869722,
               232.01846491, 232.19108753, 228.4945825 , 223.78615178,
               226.04215082, 227.52878871, 226.75049881, 224.35610404,
               223.50191971, 222.43939773, 218.78604836, 215.53746927,
               211.91537054, 207.35724055, 202.86905248, 199.08474802,
               196.16950915, 194.18881342, 192.63372173, 190.95809183,
               189.60240903, 188.46994513, 187.35236249, 186.22436297,
               185.28684358, 184.62611563, 184.5159943 , 184.96392024,
               185.43714431, 185.36720239, 184.89844269, 184.60230561,
               184.65439002, 184.85082266, 185.07255344, 185.46839496])

        # fmt: on

        # Run
        invdepth, invstrength = get_stats(temp, height)

        # Assert
        assert np.isclose(invdepth, 0.052507269999999995)
        assert np.isclose(invstrength, 18.98551173)


# Run the tests (e.g. in Spyder):
testGetStats = TestGetStats()
testGetStats.test_one_inv()
testGetStats.test_mini_noninv_near_surf()
testGetStats.test_mini_noninv_near_surf_inv_to_top
testGetStats.test_several_switches()
testGetStats.test_inv_to_top_with_thick_noninv()
testGetStats.test_simple_inv()
testGetStats.test_all_inv()
testGetStats.test_inv_goes_to_top()
testGetStats.test_one_non_inv_level()
testGetStats.test_invstrength_lt0()
