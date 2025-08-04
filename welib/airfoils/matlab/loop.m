function v = loop(v,n);
% loop ensures that v stays within 1:n
v=mod(v-1,n)+1;
