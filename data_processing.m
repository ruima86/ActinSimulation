clear;
options=2;
iii=3;  % iii runs from 1 to 5
ensemble=1;
draw = false;
if options==1
    cross_rate=1;
    L_fil=50;
    N_dendral=100;
    radius=25*(iii*0.25+0.75);
    direct='data_cylinder_radius';
elseif options==2
    cross_rate=.01;
    L_fil=iii*10;
    N_dendral=100;
    radius=25;
    direct='data_cylinder_length';
else
    cross_rate=0.01;
    L_fil=40;
    N_dendral=iii*20;
    radius=25;
    direct='data_cylinder_number';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% With Cylinder %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename1=['./',direct,'/Trajectory_',num2str(ensemble),'_crossrate_',num2str(cross_rate),...
           '_Ndendral_',num2str(N_dendral),'_Length_',num2str(L_fil),'_radius_',num2str(radius),'_date_30-Dec-2014'];
% fid1=fopen(filename1,'w');
% for tt=0:100
%     for i=1:N_dendral
%         LL=L_fil*2.7;
%         phi=(i-0.5)*2*pi/N_dendral;
%         CenterOfMass=[LL/2*cos(phi),LL/2*sin(phi),0];
%         Orientation=[sin(phi),-cos(phi),0];
%         fprintf(fid1,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f \n',tt,LL,CenterOfMass,Orientation);
%     end
% end

traj=load(filename1);

Rscope=radius+200;
angle=0:pi/100:2*pi;
xcircle=Rscope*cos(angle);
ycircle=Rscope*sin(angle);
dtheta=pi/10;
theta=-pi:dtheta:3*pi;
for t=1:1:99
    intensity=zeros(1,length(theta)-1);
    index1=find(abs(traj(:,1)-t)<1e-10);
    for j=1:length(index1)
        L=traj(index1(j),2);
        C=traj(index1(j),3:5);
        N=traj(index1(j),6:8);
        plus=C+.5*L*N;
        minus=C-.5*L*N;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Determine the start and ending positions of the filament inside the cylinder
        if norm(plus(1:2))<Rscope           
            if norm(minus(1:2))<Rscope      % plus and minus
                start=minus;ending=plus;
            else                            % plus
                % calculate the intersection point of the cylinder and the
                % line
                AA=N(1)^2+N(2)^2;
                BB=C(1)*N(1)+C(2)*N(2);
                CC=C(1)*C(1)+C(2)*C(2)-Rscope^2;
                t1=(-BB+sqrt(BB^2-AA*CC))/AA;
                t2=(-BB-sqrt(BB^2-AA*CC))/AA;
                if t1/L<0.5 && t1/L>-0.5
                    start=C+t1*N; % intersection point
                    ending=plus;
                elseif t2/L<0.5 && t2/L>-0.5
                    start=C+t2*N; % intersection point
                    ending=plus;
                end
            end
        else
            if norm(minus(1:2))<Rscope      % minus
                % calculate the intersection point of the cylinder and the
                % line
                AA=N(1)^2+N(2)^2;
                BB=C(1)*N(1)+C(2)*N(2);
                CC=C(1)*C(1)+C(2)*C(2)-Rscope^2;
                t1=(-BB+sqrt(BB^2-AA*CC))/AA;
                t2=(-BB-sqrt(BB^2-AA*CC))/AA;
                if t1/L<0.5 && t1/L>-0.5
                    start=minus;
                    ending=C+t1*N; % intersection point
                elseif t2/L<0.5 && t2/L>-0.5
                    start=minus;
                    ending=C+t2*N; % intersection point
                end
            else                            % none
                AA=N(1)^2+N(2)^2;
                BB=C(1)*N(1)+C(2)*N(2);
                CC=C(1)*C(1)+C(2)*C(2)-Rscope^2;
                if BB^2-AA*CC>0
                    t1=(-BB+sqrt(BB^2-AA*CC))/AA;
                    t2=(-BB-sqrt(BB^2-AA*CC))/AA;
                    start=C+t1*N;
                    ending=C+t2*N;
                else
                    continue; % no intersection point
                end
            end
        end
        
        if draw
            plot(xcircle,ycircle);
            hold on;
            plot([start(1),ending(1)],[start(2),ending(2)]);
            axis equal;
        end
        
        tt=-C(2)/L/N(2); % setting y=0
        
        if tt>-0.5 && tt<0.5 && (C(1)*N(2)-C(2)*N(1))*N(2) <0 % cross the line that -pi and pi overlaps
            if atan2(start(2),start(1))>atan2(ending(2),ending(1))
                theta_min=atan2(start(2),start(1));
                theta_max=atan2(ending(2),ending(1))+2*pi;
                % start has larger angle
            else
                theta_min=atan2(ending(2),ending(1));
                theta_max=atan2(start(2),start(1))+2*pi;
                temp=start;start=ending;ending=temp;
            end
            
            index1_a=find(theta<theta_min,1,'last');
            index2_a=find(theta>theta_max,1,'first');
            for i=index1_a:index2_a-1
                if draw
                    line([0,Rscope*cos(theta(i+1))],[0,Rscope*sin(theta(i+1))]);
                end
                if theta(i+1)<theta_max
                    % calculate the intersection of the line and the plane
                    tt=(C(2)-tan(theta(i+1))*C(1))/(tan(theta(i+1))*N(1)-N(2))/L;
                    intsect=C+tt*L*N;
                    if draw
                        plot(start(1),start(2),'go');
                        plot(intsect(1)+0.05,intsect(2)+0.05,'ro');
                    end
                    intensity(i)=intensity(i)+norm(intsect-start);
                    start=intsect;
                else
                    if draw
                        plot(start(1),start(2),'go');
                        plot(ending(1)+0.05,ending(2)+0.05,'ro');
                    end
                    intensity(i)=intensity(i)+norm(ending-start);
                end
            end
        else
            if atan2(start(2),start(1))>atan2(ending(2),ending(1))
                theta_min=atan2(ending(2),ending(1));
                theta_max=atan2(start(2),start(1));
                temp=start;start=ending;ending=temp;
                % start has smaller angle
            else
                theta_min=atan2(start(2),start(1));
                theta_max=atan2(ending(2),ending(1));
            end
            
            index1_a=find(theta<theta_min,1,'last');
            index2_a=find(theta>theta_max,1,'first');
            for i=index1_a:index2_a-1
                if draw
                    line([0,Rscope*cos(theta(i+1))],[0,Rscope*sin(theta(i+1))]);
                end
                if theta(i+1)<theta_max
                    % calculate the intersection of the line and the plane
                    tt=(C(2)-tan(theta(i+1))*C(1))/(tan(theta(i+1))*N(1)-N(2))/L;
                    intsect=C+tt*L*N;
                    if draw
                        plot(start(1),start(2),'go');
                        plot(intsect(1)+0.05,intsect(2)+0.05,'ro');
                    end
                    intensity(i)=intensity(i)+norm(intsect-start);
                    start=intsect;
                else
                    if draw
                        plot(start(1),start(2),'go');
                        plot(ending(1)+0.05,ending(2)+0.05,'ro');
                    end
                    intensity(i)=intensity(i)+norm(ending-start);
                end
            end
        end
        hold off
    end
    Len=length(intensity);
    intensity(1:Len/2)=intensity(1:Len/2)+intensity(Len/2+1:Len);
    plot(-pi+dtheta/2:dtheta:pi-dtheta/2,intensity(1:Len/2)/max(intensity(1:Len/2)));
    ylim([0,1.5])
    xlim([-pi,pi])
    set(gca,'xtick',-pi:pi/2:pi);
    set(gca,'xticklabel',{'-pi','-pi/2','0','pi/2','pi'});
    M(t+1)=getframe;
    %ylim([0,200])
end
%movie2avi(M, 'ActinStructure.avi')
