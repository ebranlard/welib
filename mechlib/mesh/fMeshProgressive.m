function x=fMeshProgressive(x_range,x_high,h_high,p)

x_high=x_high(1):h_high:x_high(2);

x_sup=fMeshProgressiveExtent(x_high(end), x_range(2), h_high, p);
x_inf=fMeshProgressiveExtent(x_high(1),   x_range(1), h_high, p);

x=[x_inf x_high x_sup];
