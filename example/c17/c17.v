module c17 (
// modify {nx1, nx7, nx3} -> [2:0] test
test[2],
test[1],
test[0],
nx2,
nx6,
nx23,
nx22);

// Start PIs
// input [2:0] test;
// input nx2;
// input nx6;
input \test[2], \test[1], \test[0], nx2, nx6;

// Start POs
// output nx23;
// output nx22;
output nx23, nx22;

// Start wires
// wire net_1;
// wire nx23;
// wire [2:0] test;
// wire net_2;
// wire nx22;
// wire nx6;
// wire net_0;
// wire net_3;
// wire nx2;
wire net_1, nx23, net_2, nx22, nx6, net_0, net_3, nx2;
wire [2:0] \test;



// Start cells
NAND2_X1 inst_5 ( .A2(net_3), .A1(net_0), .ZN(nx22) );
NAND2_X1 inst_2 ( .ZN(net_2), .A2(net_1), .A1(\test [1]) );
NAND2_X1 inst_1 ( .ZN(net_0), .A2(\test [0]), .A1(\test [2]) );
NAND2_X1 inst_4 ( .A1(net_3), .A2(net_2), .ZN(nx23) );
NAND2_X1 inst_3 ( .ZN(net_3), .A2(net_1), .A1(nx2) );
NAND2_X1 inst_0 ( .ZN(net_1), .A2(nx6), .A1(\test [0]) );

endmodule
