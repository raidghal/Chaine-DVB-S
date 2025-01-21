trellis = poly2trellis(3, [5,7]);
msg = [0 1 0 1 1 0 1];
code = convenc(msg,trellis);
tb=2;
decoded = vitdec(code,trellis,tb,'trunc','hard');


