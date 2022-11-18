data = """Infinity -> 0
0 -> 93
1 -> 103
2 -> 9
3 -> 88
4 -> 86
5 -> 19
6 -> 77
7 -> 121
8 -> 37
9 -> 56
10 -> 45
11 -> 23
12 -> 27
13 -> 2
14 -> 34
15 -> 68
16 -> 24
17 -> 89
18 -> 43
19 -> 99
20 -> 58
21 -> 16
22 -> 4
23 -> 59
24 -> 66
25 -> 95
26 -> 65
27 -> 79
28 -> 104
29 -> 28
30 -> 13
31 -> 62
32 -> 26
33 -> 87
34 -> 1
35 -> 109
36 -> 6
37 -> 81
38 -> 116
39 -> 112
40 -> 61
41 -> 63
42 -> 90
43 -> 110
44 -> 40
45 -> 32
46 -> 5
47 -> 97
48 -> 105
49 -> 17
50 -> 101
51 -> 49
52 -> 55
53 -> 117
54 -> 100
55 -> 115
56 -> 30
57 -> 41
58 -> 18
59 -> 69
60 -> 11
61 -> 33
62 -> Infinity
63 -> 96
64 -> 75
65 -> 10
66 -> 84
67 -> 108
68 -> 98
69 -> 60
70 -> 46
71 -> 64
72 -> 3
73 -> 122
74 -> 51
75 -> 92
76 -> 57
77 -> 50
78 -> 83
79 -> 111
80 -> 120
81 -> 67
82 -> 48
83 -> 22
84 -> 21
85 -> 73
86 -> 78
87 -> 44
88 -> 94
89 -> 74
90 -> 91
91 -> 54
92 -> 118
93 -> 31
94 -> 107
95 -> 123
96 -> 76
97 -> 52
98 -> 39
99 -> 70
100 -> 42
101 -> 36
102 -> 106
103 -> 119
104 -> 38
105 -> 80
106 -> 25
107 -> 72
108 -> 8
109 -> 53
110 -> 20
111 -> 113
112 -> 15
113 -> 12
114 -> 35
115 -> 47
116 -> 29
117 -> 114
118 -> 71
119 -> 14
120 -> 82
121 -> 85
122 -> 7
123 -> 102"""

rows = data.split("\n")
frows = [rows[i::17] for i in range(17)]

delrows = [" & ".join(row).replace("->", "&") for row in frows]

print(" \\\\\n".join(delrows))