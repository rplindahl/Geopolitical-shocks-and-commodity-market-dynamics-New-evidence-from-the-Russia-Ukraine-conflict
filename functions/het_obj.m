function f = het_obj(x,S,V,K)
   y   = x(1:K);
   c_w = x(K+1);
   f   = (S - (vech(y*y')*c_w))'/V*(S - (vech(y*y')*c_w));
end