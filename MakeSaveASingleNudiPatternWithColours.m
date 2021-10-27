%Add colours to a single nudibranch colour pattern outline in a .mat fle
clear all;
close all;
[fname,path]=uigetfile('.mat','Select a single outline .mat file');
load(fname);
nls=max(zcodes); %get number of different codes for each pattern

%VARIABLES **********************************************************************
% fname      pattern name with .mat
% path       where the file is
% outline{n} patch outlines to make zone maps
% zcode      zone codes as order of laying out the zone map (smallest later)
% nls        number of different codes/color patch classes
%********************************************************************************
%show the design  
name=fname; edges=outlines; codes=zcodes; np=length(codes);
figure; hold on;
for p=1:np
  edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
  plot(x,y,'k'); 
  text(x(1),y(1),num2str(cd)); %ring label
end;
axis equal; axis image; axis off; 

%make it
%dialog box for entering RGB then do it
ttl='ENTER RGB VALUES FOR EACH LAYER AS R,G,B';
switch nls
  case 1
    text(220,450,'1','FontWeight','bold','FontSize',14);    
    text(169,376,'2','FontWeight','bold','FontSize',14); 
    prompt=cell(nls,1); prompt{1}='Layer 1 R,G,B  SEPARATED by ,';
    for k=2:nls prompt{k}=['Layer ' num2str(k) ' R,G,B  SEPARATED by ,']; end;
    clrstr=inputdlg(prompt,ttl);
    clrs=zeros(nls,3);
    for k=1:nls
      str=clrstr{k}; str=split(str,','); str=char(str); rgb=str2num(str)';
      clrs(k,:)=rgb;
    end;
    for p=1:np
      edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
      patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
    end; 
    title('');
  case 2
    text(220,450,'1','FontWeight','bold','FontSize',14);    
    text(169,390,'2','FontWeight','bold','FontSize',14); 
    prompt=cell(nls,1); prompt{1}='Layer 1 R,G,B SEPARATED BY ,';
    for k=2:nls prompt{k}=['Layer ' num2str(k) ' R,G,B SEPARATED BY ,']; end;
    clrstr=inputdlg(prompt,ttl);
    clrs=zeros(nls,3);
    for k=1:nls
      str=clrstr{k}; str=split(str,','); str=char(str); rgb=str2num(str)';
      clrs(k,:)=rgb;
    end;
    for p=1:np
      edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
      patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
    end;  
    title('');
  case 3
    text(150,425,'1','FontWeight','bold','FontSize',14); 
    text(205,280,'2','FontWeight','bold','FontSize',14);
    text(230,324,'3','FontWeight','bold','FontSize',14);
     prompt=cell(nls,1); prompt{1}='Layer 1 R,G,B SEPARATED BY ,';
    for k=2:nls prompt{k}=['Layer ' num2str(k) ' R,G,B SEPARATED BY ,']; end;
    clrstr=inputdlg(prompt,ttl);
    clrs=zeros(nls,3);
    for k=1:nls
      str=clrstr{k}; str=split(str,','); str=char(str); rgb=str2num(str)';
      clrs(k,:)=rgb;
    end;
    for p=1:np
      edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
      patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
    end;  
    title('');  
  case 4
    text(335,425,'1','FontWeight','bold','FontSize',14); 
    text(375,386,'2','FontWeight','bold','FontSize',14);
    text(455,392,'3','FontWeight','bold','FontSize',14);
    prompt=cell(nls,1); prompt{1}='Layer 1 R,G,B SEPARATED BY ,';
    for k=2:nls prompt{k}=['Layer ' num2str(k) ' R,G,B SEPARATED BY ,']; end;
    clrstr=inputdlg(prompt,ttl);
    clrs=zeros(nls,3);
    for k=1:nls
      str=clrstr{k}; str=split(str,','); str=char(str); rgb=str2num(str)';
      clrs(k,:)=rgb;
    end;
    for p=1:np
      edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
      patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
    end;  
    title(''); 
  case 5  
    text(275,325,'1','FontWeight','bold','FontSize',14); 
    text(325,325,'2','FontWeight','bold','FontSize',14);
    text(410,340,'3','FontWeight','bold','FontSize',14);
    plot([393 408],[330 338],'b','LineWidth',3);
    text(364,327,'4','FontWeight','bold','FontSize',14);
    prompt=cell(nls,1); prompt{1}='Layer 1 R,G,B SEPARATED BY ,';
    for k=2:nls prompt{k}=['Layer ' num2str(k) ' R,G,B SEPARATED BY ,']; end;
    clrstr=inputdlg(prompt,ttl);
    clrs=zeros(nls,3);
    for k=1:nls
      str=clrstr{k}; str=split(str,','); str=char(str); rgb=str2num(str)';
      clrs(k,:)=rgb;
    end;
    for p=1:np
      edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
      patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
    end;  
    title(''); 
  case 6   
    text(400,486,'1','FontWeight','bold','FontSize',14); 
    text(200,410,'2','FontWeight','bold','FontSize',14);
    text(760,375,'2','FontWeight','bold','FontSize',14);
    text(300,445,'3','FontWeight','bold','FontSize',14);    
    text(810,340,'3','FontWeight','bold','FontSize',14);
    prompt=cell(nls,1); prompt{1}='Layer 1 RIM: R,G,B SEPARATED BY ,';
    for k=2:nls prompt{k}=['Layer ' num2str(k) ' R,G,B SEPARATED BY ,']; end;
    clrstr=inputdlg(prompt,ttl);
    clrs=zeros(nls,3);
    for k=1:nls
      str=clrstr{k}; str=split(str,','); str=char(str); rgb=str2num(str)';
      clrs(k,:)=rgb;
    end;
    for p=1:np
      edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
      patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
    end;  
    title(''); 
  case 7
    text(400,492,'1','FontWeight','bold','FontSize',14); 
    text(450,440,'2','FontWeight','bold','FontSize',14);
    text(800,420,'2','FontWeight','bold','FontSize',14);
    text(350,358,'3','FontWeight','bold','FontSize',14); 
    text(810,358,'3','FontWeight','bold','FontSize',14); 
    text(220,355,'2','FontWeight','bold','FontSize',14); 
    prompt=cell((nls-1),1); prompt{1}='Layer 1 RIM: R,G,B SEPARATED BY ,';
    for k=2:3 prompt{k}=['Layer ' num2str(k) ' R,G,B SEPARATED BY ,']; end;
    clrstr=inputdlg(prompt,ttl);
    clrs=zeros(nls,3);
    for k=1:(nls-1)
      str=clrstr{k}; str=split(str,','); str=char(str); rgb=str2num(str)';
      clrs(k,:)=rgb;
    end;
    clrs(4,:)=clrs(2,:);
    for p=1:np
      edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
      patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
    end;  
    title(''); 
 case 8 
    text(375,492,'1','FontWeight','bold','FontSize',14); 
    text(375,460,'2','FontWeight','bold','FontSize',14);
    text(760,350,'2','FontWeight','bold','FontSize',14);
    text(375,420,'3','FontWeight','bold','FontSize',14);    
    text(730,408,'3','FontWeight','bold','FontSize',14);
    prompt=cell(nls,1); prompt{1}='Layer 1 RIM: R,G,B SEPARATED BY ,';
    for k=2:nls prompt{k}=['Layer ' num2str(k) ' R,G,B SEPARATED BY ,']; end;
    clrstr=inputdlg(prompt,ttl);
    clrs=zeros(nls,3);
    for k=1:nls
      str=clrstr{k}; str=split(str,','); str=char(str); rgb=str2num(str)';
      clrs(k,:)=rgb;
    end;
    for p=1:np
      edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
      patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
    end;  
    title('');    
  case 9 
    text(450,506,'1','FontWeight','bold','FontSize',14); 
    text(450,475,'2','FontWeight','bold','FontSize',14);
    text(450,445,'1','FontWeight','bold','FontSize',14);
    text(450,410,'2','FontWeight','bold','FontSize',14);    
    text(450,362,'1','FontWeight','bold','FontSize',14);
    prompt={'Layer 1 RIM: R,G,B SEPARATED BY ,' 'Layer 2 R,G,B SEPARATED BY ,'};
    clrstr=inputdlg(prompt,ttl);
    clrs=zeros(nls,3);
    for k=1:2
      str=clrstr{k}; str=split(str,','); str=char(str); rgb=str2num(str)';
      clrs(k,:)=rgb;
    end;
    clrs(3,:)=clrs(1,:);
    clrs(4,:)=clrs(2,:);
    clrs(5,:)=clrs(1,:);
    for p=1:np
      edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
      patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
    end;  
    title('');     
 case 10 
    h1=text(400,515,'1','FontWeight','bold','FontSize',14); 
    h2=text(400,487,'2','FontWeight','bold','FontSize',14);
    h3=text(880,440,'2','FontWeight','bold','FontSize',14);
    h4=plot([870 878],[420 438],'b','LineWidth',3);
    h5=text(700,455,'3','FontWeight','bold','FontSize',14);
    h6=plot([682 700],[430 450],'b','LineWidth',3);
    h7=text(850,460,'3','FontWeight','bold','FontSize',14);
    h8=plot([832 849],[440 458],'b','LineWidth',3);
    prompt={'Layer 1 RIM: R,G,B SEPARATED BY ,' 'Layer 2 R,G,B SEPARATED BY ,'...
            'Layer 3 R,G,B SEPARATED BY ,'};
    clrstr=inputdlg(prompt,ttl);
    clrs=zeros(nls,3);
    for k=1:3
      str=clrstr{k}; str=split(str,','); str=char(str); rgb=str2num(str)';
      clrs(k,:)=rgb;
    end;
    clrs(4,:)=clrs(2,:);
    clrs(5,:)=clrs(3,:);
    delete(h1); delete(h2); delete(h3); delete(h4); delete(h5);
    delete(h6); delete(h7); delete(h8);
    for p=1:np
      edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
      patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
    end;  
    title('');     
end;
clear ax cd ay clrstr definput edg f files idx k ncds nf nls p prompt rgb str t*
clear x y zc zcodes outlines CONTENTS;

%save the outlines and colors  edges, codes, clrs
prompt={'Enter stimulus file name plus unique ID'}; ttl='STIMULUS FILE NAME';
definput={name};
fn=inputdlg(prompt,ttl,[1 25],definput); fn=char(fn); fn=['STIM' fn];
clear prompt ttl definput;
CONTENTS={'name       name of design' ...
    'np         number of patch outlines' ...
    'edges{np}  patch outline x,y coordinates' ...
    'codes(np)  patch colour codes' ...
    'clrs(nc,3) RGB values for each color'}; CONTENTS=CONTENTS';
save(fn);
fprintf(1,'Saved design in %s.mat\n',fn);

clear;