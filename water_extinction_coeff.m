function water = water_extinction_coeff

% Taken from http://www.ucl.ac.uk/medphys/research/borl/intro/spectra
% Extinction spectrum measured by Matcher 

%The extinction coefficient of water at 37 degrees over 600 - 1050 nm. The
%extinction coefficient k relates to the absorption coefficient as mu_a =
%ln(10)*k
%extinction coefficient in OD/cm.
water=[599.000	0.00076
600.000	0.00088
601.000	0.00099
602.000	0.00106
603.000	0.00105
604.000	0.00105
605.000	0.00108
606.000	0.00108
607.000	0.00108
608.000	0.00112
609.000	0.00110
610.000	0.00117
611.000	0.00121
612.000	0.00122
613.000	0.00126
614.000	0.00117
615.000	0.00128
616.000	0.00122
617.000	0.00135
618.000	0.00133
619.000	0.00131
620.000	0.00126
621.000	0.00129
622.000	0.00128
623.000	0.00135
624.000	0.00127
625.000	0.00133
626.000	0.00126
627.000	0.00125
628.000	0.00124
629.000	0.00137
630.000	0.00132
631.000	0.00135
632.000	0.00118
633.000	0.00126
634.000	0.00127
635.000	0.00114
636.000	0.00123
637.000	0.00119
638.000	0.00130
639.000	0.00130
640.000	0.00137
641.000	0.00124
642.000	0.00121
643.000	0.00140
644.000	0.00129
645.000	0.00137
646.000	0.00130
647.000	0.00139
648.000	0.00139
649.000	0.00136
650.000	0.00139
651.000	0.00140
652.000	0.00145
653.000	0.00144
654.000	0.00147
655.000	0.00141
656.000	0.00144
657.000	0.00142
658.000	0.00159
659.000	0.00153
660.000	0.00166
661.000	0.00166
662.000	0.00163
663.000	0.00164
664.000	0.00159
665.000	0.00164
666.000	0.00164
667.000	0.00164
668.000	0.00161
669.000	0.00167
670.000	0.00166
671.000	0.00162
672.000	0.00168
673.000	0.00171
674.000	0.00178
675.000	0.00166
676.000	0.00166
677.000	0.00171
678.000	0.00172
679.000	0.00172
680.000	0.00177
681.000	0.00183
682.000	0.00190
683.000	0.00183
684.000	0.00181
685.000	0.00188
686.000	0.00193
687.000	0.00200
688.000	0.00204
689.000	0.00212
690.000	0.00205
691.000	0.00220
692.000	0.00225
693.000	0.00224
694.000	0.00228
695.000	0.00238
696.000	0.00252
697.000	0.00257
698.000	0.00260
699.000	0.00278
700.000	0.00285
701.000	0.00280
702.000	0.00283
703.000	0.00296
704.000	0.00301
705.000	0.00305
706.000	0.00315
707.000	0.00338
708.000	0.00344
709.000	0.00347
710.000	0.00365
711.000	0.00388
712.000	0.00395
713.000	0.00415
714.000	0.00435
715.000	0.00457
716.000	0.00475
717.000	0.00494
718.000	0.00522
719.000	0.00552
720.000	0.00561
721.000	0.00592
722.000	0.00622
723.000	0.00651
724.000	0.00666
725.000	0.00707
726.000	0.00739
727.000	0.00782
728.000	0.00820
729.000	0.00861
730.000	0.00911
731.000	0.00969
732.000	0.01025
733.000	0.01087
734.000	0.01137
735.000	0.01182
736.000	0.01220
737.000	0.01250
738.000	0.01274
739.000	0.01290
740.000	0.01314
741.000	0.01310
742.000	0.01315
743.000	0.01320
744.000	0.01325
745.000	0.01312
746.000	0.01311
747.000	0.01308
748.000	0.01304
749.000	0.01312
750.000	0.01296
751.000	0.01306
752.000	0.01306
753.000	0.01307
754.000	0.01304
755.000	0.01297
756.000	0.01288
757.000	0.01296
758.000	0.01289
759.000	0.01278
760.000	0.01275
761.000	0.01277
762.000	0.01261
763.000	0.01252
764.000	0.01251
765.000	0.01245
766.000	0.01249
767.000	0.01239
768.000	0.01227
769.000	0.01217
770.000	0.01218
771.000	0.01208
772.000	0.01201
773.000	0.01185
774.000	0.01183
775.000	0.01191
776.000	0.01170
777.000	0.01154
778.000	0.01149
779.000	0.01150
780.000	0.01142
781.000	0.01127
782.000	0.01118
783.000	0.01117
784.000	0.01099
785.000	0.01087
786.000	0.01072
787.000	0.01055
788.000	0.01061
789.000	0.01053
790.000	0.01036
791.000	0.01024
792.000	0.01022
793.000	0.01012
794.000	0.01005
795.000	0.00983
796.000	0.00974
797.000	0.00969
798.000	0.00964
799.000	0.00951
800.000	0.00956
801.000	0.00951
802.000	0.00936
803.000	0.00939
804.000	0.00935
805.000	0.00932
806.000	0.00928
807.000	0.00923
808.000	0.00929
809.000	0.00938
810.000	0.00942
811.000	0.00952
812.000	0.00958
813.000	0.00969
814.000	0.00964
815.000	0.00984
816.000	0.00995
817.000	0.01011
818.000	0.01021
819.000	0.01019
820.000	0.01047
821.000	0.01056
822.000	0.01060
823.000	0.01093
824.000	0.01123
825.000	0.01144
826.000	0.01218
827.000	0.01260
828.000	0.01333
829.000	0.01384
830.000	0.01459
831.000	0.01510
832.000	0.01586
833.000	0.01638
834.000	0.01656
835.000	0.01699
836.000	0.01740
837.000	0.01740
838.000	0.01758
839.000	0.01773
840.000	0.01795
841.000	0.01813
842.000	0.01805
843.000	0.01814
844.000	0.01822
845.000	0.01835
846.000	0.01848
847.000	0.01856
848.000	0.01883
849.000	0.01898
850.000	0.01913
851.000	0.01908
852.000	0.01936
853.000	0.01931
854.000	0.01930
855.000	0.01934
856.000	0.01959
857.000	0.01960
858.000	0.01969
859.000	0.01972
860.000	0.01981
861.000	0.02001
862.000	0.02002
863.000	0.02008
864.000	0.01996
865.000	0.02034
866.000	0.02021
867.000	0.02045
868.000	0.02044
869.000	0.02085
870.000	0.02108
871.000	0.02116
872.000	0.02111
873.000	0.02138
874.000	0.02160
875.000	0.02172
876.000	0.02174
877.000	0.02217
878.000	0.02236
879.000	0.02278
880.000	0.02282
881.000	0.02290
882.000	0.02355
883.000	0.02366
884.000	0.02413
885.000	0.02421
886.000	0.02437
887.000	0.02440
888.000	0.02447
889.000	0.02504
890.000	0.02592
891.000	0.02579
892.000	0.02573
893.000	0.02612
894.000	0.02602
895.000	0.02643
896.000	0.02630
897.000	0.02667
898.000	0.02740
899.000	0.02781
900.000	0.02858
901.000	0.02869
902.000	0.02909
903.000	0.02938
904.000	0.02962
905.000	0.03003
906.000	0.03042
907.000	0.03073
908.000	0.03106
909.000	0.03155
910.000	0.03200
911.000	0.03236
912.000	0.03294
913.000	0.03352
914.000	0.03420
915.000	0.03476
916.000	0.03558
917.000	0.03642
918.000	0.03731
919.000	0.03822
920.000	0.03941
921.000	0.04067
922.000	0.04191
923.000	0.04338
924.000	0.04499
925.000	0.04674
926.000	0.04854
927.000	0.05052
928.000	0.05266
929.000	0.05490
930.000	0.05740
931.000	0.05980
932.000	0.06248
933.000	0.06530
934.000	0.06820
935.000	0.07118
936.000	0.07426
937.000	0.07749
938.000	0.08087
939.000	0.08434
940.000	0.08804
941.000	0.09178
942.000	0.09607
943.000	0.10042
944.000	0.10503
945.000	0.11001
946.000	0.11512
947.000	0.12066
948.000	0.12677
949.000	0.13335
950.000	0.14015
951.000	0.14739
952.000	0.15492
953.000	0.16251
954.000	0.17015
955.000	0.17750
956.000	0.18438
957.000	0.19066
958.000	0.19641
959.000	0.20144
960.000	0.20562
961.000	0.20918
962.000	0.21194
963.000	0.21397
964.000	0.21561
965.000	0.21687
966.000	0.21790
967.000	0.21857
968.000	0.21913
969.000	0.21948
970.000	0.21969
971.000	0.21980
972.000	0.21992
973.000	0.21976
974.000	0.21959
975.000	0.21926
976.000	0.21863
977.000	0.21796
978.000	0.21714
979.000	0.21609
980.000	0.21491
981.000	0.21361
982.000	0.21211
983.000	0.21059
984.000	0.20892
985.000	0.20710
986.000	0.20524
987.000	0.20322
988.000	0.20111
989.000	0.19895
990.000	0.19675
991.000	0.19451
992.000	0.19217
993.000	0.18985
994.000	0.18732
995.000	0.18488
996.000	0.18242
997.000	0.17985
998.000	0.17735
999.000	0.17484
1000.000	0.17227
1001.000	0.16973
1002.000	0.16703
1003.000	0.16442
1004.000	0.16173
1005.000	0.15907
1006.000	0.15656
1007.000	0.15394
1008.000	0.15120
1009.000	0.14863
1010.000	0.14602
1011.000	0.14320
1012.000	0.14060
1013.000	0.13796
1014.000	0.13537
1015.000	0.13282
1016.000	0.13020
1017.000	0.12760
1018.000	0.12510
1019.000	0.12270
1020.000	0.12018
1021.000	0.11773
1022.000	0.11526
1023.000	0.11294
1024.000	0.11059
1025.000	0.10822
1026.000	0.10583
1027.000	0.10352
1028.000	0.10120
1029.000	0.09910
1030.000	0.09692
1031.000	0.09482
1032.000	0.09276
1033.000	0.09080
1034.000	0.08901
1035.000	0.08721
1036.000	0.08547
1037.000	0.08378
1038.000	0.08222
1039.000	0.08066
1040.000	0.07927
1041.000	0.07790
1042.000	0.07651
1043.000	0.07527
1044.000	0.07397
1045.000	0.07273
1046.000	0.07158
1047.000	0.07041
1048.000	0.06931
1049.000	0.06823
1050.000	0.06724
1051.000	0.06631
];

