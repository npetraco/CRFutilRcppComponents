

tm <- array(NaN,c(2,2))
testemptymat(tm)

testemptymat(array(-1,c(1,1)))

ta <- array(0,c(2,2,1))
tm <- array(0,c(2,2))
ta
tm
testslicedetect2(tm)
testslicedetect2(ta)

testslicedetect3(tm)
testslicedetect3(ta)

testslicedetect4(tm)
testslicedetect4(ta)

testslicedetect5(tm)
testslicedetect5(ta)
