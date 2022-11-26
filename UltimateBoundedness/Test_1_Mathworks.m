close all;
clear all;
clc

% create a fig file
     fnam='Fig1.fig';
     fh=figure;
     surfl(peaks(32));
     shading interp;
     camlight right;
     lighting phong;
     title('TEST');
     saveas(gcf,fnam);
     delete(fh);
% the engine
% ...a subplot template
     fhs=figure;
for i=1:2
     sh(i)=subplot(2,2,i);
end
% ...now, create the current fig
     axes(sh(1));
     plot(rand(10,1));
% ...and, reload the fig file
     fh=open(fnam);
% ...move it to the subplot
     ch=copyobj(gca,fhs);
% ...resize it
     set(ch,'position',get(sh(2),'position'));
% ...and delete the fig's canvas
     delete(fh);
% ...and the template
     delete(sh(2));