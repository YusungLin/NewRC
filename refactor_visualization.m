function visualParam = refactor_visualization(param, dtmax)
    % load parameters
    handles = param.handles;
    concretetype = param.concretetype;
    
    % 畫斷面圖 sectionView
    scoordinate = sectionView(param, dtmax);

    % MANDER Model
    if concretetype==0||concretetype==1
        mander = manderModel(param);
        visualParam.mander = mander;
    % HUNG Model
    elseif concretetype==2
        hungModel(param);
    end
    
    steelModel(param)
    
    handles.axes2.XTickLabel = string(handles.axes2.XTick);
    plot(handles.axes2, 0, 0, 'r');
    handles.axes3.YTickLabel = string(handles.axes3.YTick);
    plot(handles.axes3,0,0,'r');
    
    visualParam.handles = handles;
    visualParam.scoordinate = scoordinate;
end

function scoordinate = sectionView(param, dtmax)
    handles = param.handles;
    h = param.h;
    b = param.b;
    s = param.s;
    number = param.number;
    mat = param.mat;
    
    x=0:0.1:b;
    y=0:0.1:h;
    axes(handles.axes1)
%     if b>=h
%         axis([0 b 0 b]);
%     else
%         axis([0 h 0 h]);
%     end
    axis([0 max(b, h) 0 max(b,h)]);
    
    set(gca,'XColor','k');
    set(gca,'YColor','k');
%     xlabel('b','color','k','FontSize',10,'FontWeight','bold');
%     ylabel('h','color','k','FontSize',10,'FontWeight','bold');
    axis equal;
    plot(x,0,'k.');
    plot(x,h,'k.');
    plot(0,y,'k.');
    plot(b,y,'k.');
    axis equal;

%     for i=1:1:mat(1)
%         number(i)=number(i)-2;
%         if number(i)>0
%             sxfirst(i,1)=h-dtmax;
%             sx(i,1)=sxfirst(i,1);
%             syfirst(1,1)=sxfirst(1,1);
%             if number(i)> 1
%                 sxend(i,(number(i)))=b-sxfirst(i,1);
%                 sx(i,(number(i)))=sxend(i,(number(i)));
%             end
%             sy(i)=s(i,1);
%             axes(handles.axes1)
%             plot(sxfirst(i,1),h-sy(i),'ko');hold on;%%%%%%
%             if (number(i))> 1
%                 plot(sxend(i,(number(i))),h-sy(i),'ko');hold on;%%%%%%
%             end
%             if (number(i))>2
%                 for j=2:1:(number(i)-1)
%                     sx(i,j)=sxfirst(i,1)+(j-1)*(((sxend(i,(number(i))))-sxfirst(i,1))/(number(i)-1));
%                     axes(handles.axes1)
%                     plot(sx(i,j),h-sy(i),'ko');hold on;%%%%%
%                 end
%             end
%         end
%     end
    number = number - 2;
    sxfirst = h - dtmax;
    sxend = b - sxfirst;
    sy = s(:,1)';
    for i = 1:mat(1)
        if number(i) > 0
            sx(i, number(i)) = sxend;
            sx(i, 1) = sxfirst;
            for j = 2:(number(i) - 1)
                sx(i,j) = sxfirst + (j - 1) * (sxend - sxfirst) / (number(i) - 1);
                plot(sx(i, j), h - sy(i), 'ko');
            end
            plot(sx(i, number(i)), h - sy(i), 'ko');
            plot(sx(i, 1), h - sy(i), 'ko');
        end
    end

    p = 1;
    for i = 1:mat(1)
        for f = 1:(number(i))
            scoordinate(:,p) = [sx(i,f), sy(i)];
            p = p + 1;
        end
    end

%     betal=1.05-(0.05*(fc/70));
%     if betal<0.65
%         betal=0.65;
%     end
%     if betal>0.85
%         betal=0.85;
%     end
end

function mander = manderModel(param)
    % unpack parameters
    handles = param.handles;
    fc = param.mander.fc;
    Ec = param.mander.Ec;
    K = param.mander.K;
    ecu = param.mander.ecu;
    ft = param.mander.ft;

    unfcc = fc;
    if unfcc <= 204
        unecc = 0.002;
    elseif unfcc > 204 || unfcc < 1020
        unecc = 0.002 + 0.001 * (unfcc - 204) / 816;
    else
        unecc = 0.003;
    end
    unEsec = unfcc / unecc;
    unr = Ec / (Ec - unEsec);
    fcc = K * fc;
    ecc = unecc * (1 + 5 * (K - 1));
    Esec = fcc / ecc;
    r = Ec / (Ec - Esec);

%     v=1;u=1;
%     for xx=0:0.0001:ecu
%         confinedyy(u,v)= (r*(xx/ecc))/((r-1)+(xx/ecc).^r)*fcc;
%         v=v+1;
%     end
%     v=1;u=1;
%     for xx=0:0.0001:0.004
%         unconfinedyy(u,v)= (unr*(xx/unecc))/((unr-1)+(xx/unecc).^unr)*unfcc;
%         v=v+1;
%     end
%     unconfined=unconfinedyy(u,v-1);
    x1 = 0:0.0001:ecu;
    confinedyy = (r * (x1 / ecc)) ./ ((r - 1) + (x1 / ecc) .^ r) * fcc;
    x2 = 0:0.0001:0.004;
    unconfinedyy = (unr * (x2 / unecc)) ./ ((unr - 1) + (x2 / unecc) .^ unr) * unfcc;
    unconfined = unconfinedyy(end);
    ectension = ft / Ec;
    
    axes(handles.axes2);
    xtickformat('%.3f');
    grid on;
    
    plot(x1, confinedyy, 'r');
    plot(x2, unconfinedyy, 'r');
    line([0.004,0.006], [unconfined,0], 'Color', 'red');
    line([0,-ectension], [0,-ft], 'Color', 'blue');

    mander.unr = unr;
    mander.unecc = unecc;
    mander.unfcc = unfcc;
    mander.r = r;
    mander.ecc = ecc;
    mander.fcc = fcc;
    mander.ectension = ectension;
    mander.unconfined = unconfined;
end

function hungModel(param)
    % unpack parameters
    SIGtc = param.hung.SIGtc;
    SIGtp = param.hung.SIGtp;
    EPSILONtc = param.hung.EPSILONtc;
    EPSILONtp = param.hung.EPSILONtp;
    EPSILONtu = param.hung.EPSILONtu;
    SIGcp = param.hung.SIGcp;
    SIGcu = param.hung.SIGcu;
    EPSILONcp = param.hung.EPSILONcp;
    EPSILONcu = param.hung.EPSILONcu;
    
%     v=1;u=1;
%     for xx=0:0.0001:EPSILONcp
%         HPFRCCyy(u,v)= SIGcp*(2*(xx/EPSILONcp)-(xx/EPSILONcp)^2);
%         v=v+1;
%     end
    x = 0:0.0001:EPSILONcp;
    HPFRCCyy = SIGcp * (2 * (x / EPSILONcp) - (x / EPSILONcp) ^ 2);
    
    axes(handles.axes2)
    grid on;
    plot(-xx,-HPFRCCyy,'r');
    xtickformat('%.3f');
    
    line([0,EPSILONtc],[0,SIGtc],'Color','blue');
    line([EPSILONtc,EPSILONtp],[SIGtc,SIGtp],'Color','blue');
    line([EPSILONtp,EPSILONtu],[SIGtp,0],'Color','blue');
    line([-EPSILONcp,-EPSILONcu],[-SIGcp,-SIGcu],'Color','red');
    line([-EPSILONcu,-0.1],[-SIGcu,-SIGcu],'Color','red');
end

function steelModel(param)
    % unpack parameters
    handles = param.handles;
    stnum = param.stnum;
    fy = param.fy;
    fu = param.fu;
    Es = param.Es;
    esh = param.esh;
    esu = param.esu;
    power = param.power;
    fcr = param.fcr;
    
    axes(handles.axes3);
    grid on;
    
%     for num=1:1:stnum
%         steely=[];
%         axes(handles.axes3)
%         line([0,fy(num)/Es(num)],[0,fy(num)],'Color','blue');hold on;
%         line([fy(num)/Es(num),esh(num)],[fy(num),fy(num)],'Color','blue');hold on;
%         f=1;
%         for x=esh(num):0.01:esu(num)
%             steely(1,f)=fu(num)+(fy(num)-fu(num))*((esu(num)-x)/(esu(num)-esh(num))).^power(num);
%             f=f+1;
%         end
%         axes(handles.axes3)
%         x=esh(num):0.01:esu(num);
%         plot(x,steely,'b');hold on;grid on;
    for num = 1:stnum
        line([0,fy(num) / Es(num)],[0,fy(num)],'Color','blue');
        line([fy(num) / Es(num),esh(num)],[fy(num),fy(num)],'Color','blue');
        
        x = esh(num):0.01:esu(num);
        steely = fu(num) + (fy(num) - fu(num)) * ((esu(num) - x) / (esu(num) - esh(num))) .^ power(num);
        plot(x,steely,'b');
%     end

%     for num=1:1:stnum
%         steelxx=[];steelyy=[];
%         axes(handles.axes3)
        line([0, -(fy(num) / Es(num))], [0, -fy(num)], 'Color', 'red');
        %%line([-0.004,-(fy(num)/Es(num))],[-fy(num),-fy(num)],'Color','red');hold on;
        if fcr(num) > fy(num)
            line([-(fy(num) / Es(num)), -esh(num)], [-fy(num), -fy(num)], 'Color', 'red');
        end
%         f=1;
        if fcr(num)>fu(num)
%             for x=esh(num):0.01:esu(num)
%                 steelyy(1,f)=-(fu(num)+(fy(num)-fu(num))*((esu(num)-x)/(esu(num)-esh(num))).^power(num));
%                 f=f+1;
%             end
%             axes(handles.axes3)
%             x=esh(num):0.01:esu(num);
%             plot(-x,steelyy,'r');hold on;grid on;
            x = esh(num):0.01:esu(num);
            steelyy = -(fu(num) + (fy(num) - fu(num)) * ((esu(num) - x) / (esu(num) - esh(num))) .^ power(num));
            plot(-x,steelyy,'r');

        elseif fu(num) >= fcr(num) && fcr(num) >= fy(num)
            ex = esu(num) - (((fcr(num) - fu(num)) / (fy(num) - fu(num))) ^ (1 / power(num))) * (esu(num) - esh(num));
%             for x=esh(num):0.01:ex
%                 steelxx(1,f)=-x;
%                 steelyy(1,f)=-(fu(num)+(fy(num)-fu(num))*((esu(num)-x)/(esu(num)-esh(num))).^power(num));
%                 f=f+1;
%             end
%             axes(handles.axes3)
            x = esh(num):0.01:ex;
            steelxx = -x;
            steelyy = -(fu(num) + (fy(num) - fu(num)) * ((esu(num) - x) / (esu(num) - esh(num))) .^ power(num));
            plot(steelxx,steelyy,'r');
        end
    end
end