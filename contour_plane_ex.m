%% contour plane example

[X,Y,Z] = peaks(30);
h1=surf(X,Y,Z);
xmin=-3;xmax=3;ymin=-3;ymax=3;zmin=-10;zmax=10;
axis([xmin xmax ymin ymax zmin zmax]) 
hold on;
%% contour on yz plane
[C,h]=contour(Y,Z,X);
hpatch = get(h,'Children');
for i = 1:length(hpatch)
      ch = hpatch(i);
      xdata = get(ch,'Xdata'); %
      ydata = get(ch,'Ydata');
      set(ch,'Xdata',zeros(size(xdata))+xmax);
      set(ch,'Zdata',ydata);
      set(ch,'Ydata',xdata);
end
%% contour on xz plane
[C,h]=contour(X,Z,Y);
hpatch = get(h,'Children');
for i = 1:length(hpatch)
      ch = hpatch(i);
      xdata = get(ch,'Xdata'); %
      ydata = get(ch,'Ydata');
      set(ch,'Ydata',zeros(size(ydata))+ymax);
      set(ch,'Zdata',ydata);
      set(ch,'Xdata',xdata);
end