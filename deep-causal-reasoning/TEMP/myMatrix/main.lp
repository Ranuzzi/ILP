\ENCODING=ISO-8859-1
\Problem name: 

Minimize
 obj: 100 abs_{a,0}#20 + 100 abs_{c,0}#21 + 100 abs_{b,0}#22 + 100 abs_{e,0}#23
      + 100 abs_{d,0}#24 + 100 abs_{g,0}#25 + 100 abs_{f,0}#26
      + 100 abs_{h,0}#27 + 100 abs_{j,0}#28 + 100 abs_{y,0}#29 + 20 xp_{a,0}#30
      + 20 xp_{c,0}#31 + 20 xp_{b,0}#32 + 20 xp_{e,0}#33 + 20 xp_{d,0}#34
      + 20 xp_{g,0}#35 + 20 xp_{f,0}#36 + 20 xp_{h,0}#37 + 20 xp_{j,0}#38
      + 20 xp_{y,0}#39 + 20 xm_{a,0}#40 + 20 xm_{c,0}#41 + 20 xm_{b,0}#42
      + 20 xm_{e,0}#43 + 20 xm_{d,0}#44 + 20 xm_{g,0}#45 + 20 xm_{f,0}#46
      + 20 xm_{h,0}#47 + 20 xm_{j,0}#48 + 20 xm_{y,0}#49 + 0 dist_{j,0}#82
      + 0 dist_{y,0}#83
Subject To
 c1:   - x_{a,0}#0 + up_{0,0}#50 >= 0
 c2:   - x_{a,0}#0 + up_{1,0}#51 >= 0
 c3:   - x_{b,0}#2 + up_{2,0}#52 >= 0
 c4:   - x_{c,0}#1 + up_{3,0}#53 >= 0
 c5:   - x_{d,0}#4 + up_{4,0}#54 >= 0
 c6:   - x_{e,0}#3 + up_{5,0}#55 >= 0
 c7:   x_{b,0}#2 + up_{6,0}#56 >= 0
 c8:   - x_{f,0}#6 + up_{7,0}#57 >= 0
 c9:   x_{h,0}#7 + up_{8,0}#58 >= 0
 c10:  - x_{g,0}#5 + up_{9,0}#59 >= 0
 c11:  - x_{j,0}#8 + up_{10,0}#60 >= 0
 c12:  x_{y,0}#9 + up_{11,0}#61 >= 0
 c13:  x_{a,0}#0 + um_{0,0}#62 >= 0
 c14:  x_{a,0}#0 + um_{1,0}#63 >= 0
 c15:  x_{b,0}#2 + um_{2,0}#64 >= 0
 c16:  x_{c,0}#1 + um_{3,0}#65 >= 0
 c17:  x_{d,0}#4 + um_{4,0}#66 >= 0
 c18:  x_{e,0}#3 + um_{5,0}#67 >= 0
 c19:  - x_{b,0}#2 + um_{6,0}#68 >= 0
 c20:  x_{f,0}#6 + um_{7,0}#69 >= 0
 c21:  - x_{h,0}#7 + um_{8,0}#70 >= 0
 c22:  x_{g,0}#5 + um_{9,0}#71 >= 0
 c23:  x_{j,0}#8 + um_{10,0}#72 >= 0
 c24:  - x_{y,0}#9 + um_{11,0}#73 >= 0
 c25:  up_{0,0}#50 + um_{0,0}#62 <= 1
 c26:  up_{1,0}#51 + um_{1,0}#63 <= 1
 c27:  up_{2,0}#52 + um_{2,0}#64 <= 1
 c28:  up_{3,0}#53 + um_{3,0}#65 <= 1
 c29:  up_{4,0}#54 + um_{4,0}#66 <= 1
 c30:  up_{5,0}#55 + um_{5,0}#67 <= 1
 c31:  up_{6,0}#56 + um_{6,0}#68 <= 1
 c32:  up_{7,0}#57 + um_{7,0}#69 <= 1
 c33:  up_{8,0}#58 + um_{8,0}#70 <= 1
 c34:  up_{9,0}#59 + um_{9,0}#71 <= 1
 c35:  up_{10,0}#60 + um_{10,0}#72 <= 1
 c36:  up_{11,0}#61 + um_{11,0}#73 <= 1
 c37:  x_{a,0}#0 - up_{0,0}#50 + um_{0,0}#62 <= 0
 c38:  x_{a,0}#0 - up_{1,0}#51 + um_{1,0}#63 <= 0
 c39:  x_{b,0}#2 - up_{2,0}#52 + um_{2,0}#64 <= 0
 c40:  x_{c,0}#1 - up_{3,0}#53 + um_{3,0}#65 <= 0
 c41:  x_{d,0}#4 - up_{4,0}#54 + um_{4,0}#66 <= 0
 c42:  x_{e,0}#3 - up_{5,0}#55 + um_{5,0}#67 <= 0
 c43:  - x_{b,0}#2 - up_{6,0}#56 + um_{6,0}#68 <= 0
 c44:  x_{f,0}#6 - up_{7,0}#57 + um_{7,0}#69 <= 0
 c45:  - x_{h,0}#7 - up_{8,0}#58 + um_{8,0}#70 <= 0
 c46:  x_{g,0}#5 - up_{9,0}#59 + um_{9,0}#71 <= 0
 c47:  x_{j,0}#8 - up_{10,0}#60 + um_{10,0}#72 <= 0
 c48:  - x_{y,0}#9 - up_{11,0}#61 + um_{11,0}#73 <= 0
 c49:  - x_{a,0}#0 + up_{0,0}#50 - um_{0,0}#62 <= 0
 c50:  - x_{a,0}#0 + up_{1,0}#51 - um_{1,0}#63 <= 0
 c51:  - x_{b,0}#2 + up_{2,0}#52 - um_{2,0}#64 <= 0
 c52:  - x_{c,0}#1 + up_{3,0}#53 - um_{3,0}#65 <= 0
 c53:  - x_{d,0}#4 + up_{4,0}#54 - um_{4,0}#66 <= 0
 c54:  - x_{e,0}#3 + up_{5,0}#55 - um_{5,0}#67 <= 0
 c55:  x_{b,0}#2 + up_{6,0}#56 - um_{6,0}#68 <= 0
 c56:  - x_{f,0}#6 + up_{7,0}#57 - um_{7,0}#69 <= 0
 c57:  x_{h,0}#7 + up_{8,0}#58 - um_{8,0}#70 <= 0
 c58:  - x_{g,0}#5 + up_{9,0}#59 - um_{9,0}#71 <= 0
 c59:  - x_{j,0}#8 + up_{10,0}#60 - um_{10,0}#72 <= 0
 c60:  x_{y,0}#9 + up_{11,0}#61 - um_{11,0}#73 <= 0
 c61:  xp_{a,0}#30 - up_{5,0}#55 - up_{6,0}#56 - up_{10,0}#60 <= 0
 c62:  xp_{c,0}#31 - up_{0,0}#50 <= 0
 c63:  xp_{b,0}#32 - up_{1,0}#51 - up_{7,0}#57 <= 0
 c64:  xp_{e,0}#33 - up_{3,0}#53 - up_{4,0}#54 - up_{8,0}#58 <= 0
 c65:  xp_{d,0}#34 - up_{2,0}#52 - up_{11,0}#61 <= 0
 c66:  xp_{g,0}#35 <= 0
 c67:  xp_{f,0}#36 <= 0
 c68:  xp_{h,0}#37 - up_{9,0}#59 <= 0
 c69:  xp_{j,0}#38 <= 0
 c70:  xp_{y,0}#39 <= 0
 c71:  xm_{a,0}#40 - um_{5,0}#67 - um_{6,0}#68 - um_{10,0}#72 <= 0
 c72:  xm_{c,0}#41 - um_{0,0}#62 <= 0
 c73:  xm_{b,0}#42 - um_{1,0}#63 - um_{7,0}#69 <= 0
 c74:  xm_{e,0}#43 - um_{3,0}#65 - um_{4,0}#66 - um_{8,0}#70 <= 0
 c75:  xm_{d,0}#44 - um_{2,0}#64 - um_{11,0}#73 <= 0
 c76:  xm_{g,0}#45 <= 0
 c77:  xm_{f,0}#46 <= 0
 c78:  xm_{h,0}#47 - um_{9,0}#71 <= 0
 c79:  xm_{j,0}#48 <= 0
 c80:  xm_{y,0}#49 <= 0
 c81:  x_{a,0}#0 - B_{a,0}#10 - xp_{a,0}#30 + xm_{a,0}#40  = 0
 c82:  x_{c,0}#1 - B_{c,0}#11 - xp_{c,0}#31 + xm_{c,0}#41  = 0
 c83:  x_{b,0}#2 - B_{b,0}#12 - xp_{b,0}#32 + xm_{b,0}#42  = 0
 c84:  x_{e,0}#3 - B_{e,0}#13 - xp_{e,0}#33 + xm_{e,0}#43  = 0
 c85:  x_{d,0}#4 - B_{d,0}#14 - xp_{d,0}#34 + xm_{d,0}#44  = 0
 c86:  x_{g,0}#5 - B_{g,0}#15 - xp_{g,0}#35 + xm_{g,0}#45  = 0
 c87:  x_{f,0}#6 - B_{f,0}#16 - xp_{f,0}#36 + xm_{f,0}#46  = 0
 c88:  x_{h,0}#7 - B_{h,0}#17 - xp_{h,0}#37 + xm_{h,0}#47  = 0
 c89:  x_{j,0}#8 - B_{j,0}#18 - xp_{j,0}#38 + xm_{j,0}#48  = 0
 c90:  x_{y,0}#9 - B_{y,0}#19 - xp_{y,0}#39 + xm_{y,0}#49  = 0
 c91:  x_{c,0}#1 + abs_{c,0}#21 >= -1
 c92:  x_{b,0}#2 + abs_{b,0}#22 >= -1
 c93:  x_{e,0}#3 + abs_{e,0}#23 >= 1
 c94:  - x_{c,0}#1 + abs_{c,0}#21 >= 1
 c95:  - x_{b,0}#2 + abs_{b,0}#22 >= 1
 c96:  - x_{e,0}#3 + abs_{e,0}#23 >= -1
 c97:  - 101 up_{0,0}#50 - dist_{a,0}#74 + dist_{c,0}#75 >= -100
 c98:  - 101 up_{1,0}#51 - dist_{a,0}#74 + dist_{b,0}#76 >= -100
 c99:  - 101 up_{3,0}#53 - dist_{c,0}#75 + dist_{e,0}#77 >= -100
 c100: - 101 up_{4,0}#54 + dist_{e,0}#77 - dist_{d,0}#78 >= -100
 c101: - 101 up_{7,0}#57 + dist_{b,0}#76 - dist_{f,0}#80 >= -100
 c102: - 101 up_{8,0}#58 + dist_{e,0}#77 - dist_{h,0}#81 >= -100
 c103: - 101 up_{9,0}#59 - dist_{g,0}#79 + dist_{h,0}#81 >= -100
 c104: - 101 um_{0,0}#62 - dist_{a,0}#74 + dist_{c,0}#75 >= -100
 c105: - 101 um_{1,0}#63 - dist_{a,0}#74 + dist_{b,0}#76 >= -100
 c106: - 101 um_{3,0}#65 - dist_{c,0}#75 + dist_{e,0}#77 >= -100
 c107: - 101 um_{4,0}#66 + dist_{e,0}#77 - dist_{d,0}#78 >= -100
 c108: - 101 um_{7,0}#69 + dist_{b,0}#76 - dist_{f,0}#80 >= -100
 c109: - 101 um_{8,0}#70 + dist_{e,0}#77 - dist_{h,0}#81 >= -100
 c110: - 101 um_{9,0}#71 - dist_{g,0}#79 + dist_{h,0}#81 >= -100
 c111: x_{a,0}#0  = -1
 c112: x_{g,0}#5 - B_{g,0}#15  = 0
 c113: x_{f,0}#6 - B_{f,0}#16  = 0
 c114: x_{j,0}#8 - B_{j,0}#18  = 0
 c115: x_{y,0}#9 - B_{y,0}#19  = 0
 c116: B_{c,0}#11  = 0
 c117: B_{b,0}#12  = 0
 c118: B_{e,0}#13  = 0
 c119: B_{g,0}#15  = 0
 c120: B_{f,0}#16  = 0
 c121: B_{h,0}#17  = 0
 c122: B_{j,0}#18  = 0
 c123: B_{y,0}#19  = 0
Bounds
-1 <= x_{a,0}#0 <= 1
-1 <= x_{c,0}#1 <= 1
-1 <= x_{b,0}#2 <= 1
-1 <= x_{e,0}#3 <= 1
-1 <= x_{d,0}#4 <= 1
-1 <= x_{g,0}#5 <= 1
-1 <= x_{f,0}#6 <= 1
-1 <= x_{h,0}#7 <= 1
-1 <= x_{j,0}#8 <= 1
-1 <= x_{y,0}#9 <= 1
-1 <= B_{a,0}#10 <= 1
-1 <= B_{c,0}#11 <= 1
-1 <= B_{b,0}#12 <= 1
-1 <= B_{e,0}#13 <= 1
-1 <= B_{d,0}#14 <= 1
-1 <= B_{g,0}#15 <= 1
-1 <= B_{f,0}#16 <= 1
-1 <= B_{h,0}#17 <= 1
-1 <= B_{j,0}#18 <= 1
-1 <= B_{y,0}#19 <= 1
 0 <= abs_{a,0}#20 <= 2
 0 <= abs_{c,0}#21 <= 2
 0 <= abs_{b,0}#22 <= 2
 0 <= abs_{e,0}#23 <= 2
 0 <= abs_{d,0}#24 <= 2
 0 <= abs_{g,0}#25 <= 2
 0 <= abs_{f,0}#26 <= 2
 0 <= abs_{h,0}#27 <= 2
 0 <= abs_{j,0}#28 <= 2
 0 <= abs_{y,0}#29 <= 2
 0 <= xp_{a,0}#30 <= 1
 0 <= xp_{c,0}#31 <= 1
 0 <= xp_{b,0}#32 <= 1
 0 <= xp_{e,0}#33 <= 1
 0 <= xp_{d,0}#34 <= 1
 0 <= xp_{g,0}#35 <= 1
 0 <= xp_{f,0}#36 <= 1
 0 <= xp_{h,0}#37 <= 1
 0 <= xp_{j,0}#38 <= 1
 0 <= xp_{y,0}#39 <= 1
 0 <= xm_{a,0}#40 <= 1
 0 <= xm_{c,0}#41 <= 1
 0 <= xm_{b,0}#42 <= 1
 0 <= xm_{e,0}#43 <= 1
 0 <= xm_{d,0}#44 <= 1
 0 <= xm_{g,0}#45 <= 1
 0 <= xm_{f,0}#46 <= 1
 0 <= xm_{h,0}#47 <= 1
 0 <= xm_{j,0}#48 <= 1
 0 <= xm_{y,0}#49 <= 1
 0 <= up_{0,0}#50 <= 1
 0 <= up_{1,0}#51 <= 1
 0 <= up_{2,0}#52 <= 1
 0 <= up_{3,0}#53 <= 1
 0 <= up_{4,0}#54 <= 1
 0 <= up_{5,0}#55 <= 1
 0 <= up_{6,0}#56 <= 1
 0 <= up_{7,0}#57 <= 1
 0 <= up_{8,0}#58 <= 1
 0 <= up_{9,0}#59 <= 1
 0 <= up_{10,0}#60 <= 1
 0 <= up_{11,0}#61 <= 1
 0 <= um_{0,0}#62 <= 1
 0 <= um_{1,0}#63 <= 1
 0 <= um_{2,0}#64 <= 1
 0 <= um_{3,0}#65 <= 1
 0 <= um_{4,0}#66 <= 1
 0 <= um_{5,0}#67 <= 1
 0 <= um_{6,0}#68 <= 1
 0 <= um_{7,0}#69 <= 1
 0 <= um_{8,0}#70 <= 1
 0 <= um_{9,0}#71 <= 1
 0 <= um_{10,0}#72 <= 1
 0 <= um_{11,0}#73 <= 1
 0 <= dist_{a,0}#74 <= 100
 0 <= dist_{c,0}#75 <= 100
 0 <= dist_{b,0}#76 <= 100
 0 <= dist_{e,0}#77 <= 100
 0 <= dist_{d,0}#78 <= 100
 0 <= dist_{g,0}#79 <= 100
 0 <= dist_{f,0}#80 <= 100
 0 <= dist_{h,0}#81 <= 100
 0 <= dist_{j,0}#82 <= 100
 0 <= dist_{y,0}#83 <= 100
Binaries
 xp_{a,0}#30  xp_{c,0}#31  xp_{b,0}#32  xp_{e,0}#33  xp_{d,0}#34  xp_{g,0}#35 
 xp_{f,0}#36  xp_{h,0}#37  xp_{j,0}#38  xp_{y,0}#39  xm_{a,0}#40  xm_{c,0}#41 
 xm_{b,0}#42  xm_{e,0}#43  xm_{d,0}#44  xm_{g,0}#45  xm_{f,0}#46  xm_{h,0}#47 
 xm_{j,0}#48  xm_{y,0}#49  up_{0,0}#50  up_{1,0}#51  up_{2,0}#52  up_{3,0}#53 
 up_{4,0}#54  up_{5,0}#55  up_{6,0}#56  up_{7,0}#57  up_{8,0}#58  up_{9,0}#59 
 up_{10,0}#60  up_{11,0}#61  um_{0,0}#62  um_{1,0}#63  um_{2,0}#64 
 um_{3,0}#65  um_{4,0}#66  um_{5,0}#67  um_{6,0}#68  um_{7,0}#69  um_{8,0}#70 
 um_{9,0}#71  um_{10,0}#72  um_{11,0}#73 
Generals
 x_{a,0}#0  x_{c,0}#1  x_{b,0}#2  x_{e,0}#3  x_{d,0}#4  x_{g,0}#5  x_{f,0}#6 
 x_{h,0}#7  x_{j,0}#8  x_{y,0}#9  B_{a,0}#10  B_{c,0}#11  B_{b,0}#12 
 B_{e,0}#13  B_{d,0}#14  B_{g,0}#15  B_{f,0}#16  B_{h,0}#17  B_{j,0}#18 
 B_{y,0}#19  abs_{a,0}#20  abs_{c,0}#21  abs_{b,0}#22  abs_{e,0}#23 
 abs_{d,0}#24  abs_{g,0}#25  abs_{f,0}#26  abs_{h,0}#27  abs_{j,0}#28 
 abs_{y,0}#29 
End
