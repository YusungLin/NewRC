function writeParam = refactor_write_file(param, steelParam, visualParam, ...
                    srotation, crotation, cxy, dtmax, dtnum)
    % Write files.
    %   Output the model data into files prefixed with "Mocur_"
    
    % load data
    
    handles = param.handles;
    pathname = param.pathname;
    filename = param.filename;
    concretetype = param.concretetype;
    sstyp = param.sstyp;
        
    h = param.h;
    b = param.b;
    cover = param.cover;
    P = param.P;
    
    fy = param.fy;
    fu = param.fu;
    Es = param.Es;
    esh = param.esh;
    esu = param.esu;
    power = param.power;
    fcr = param.fcr;
    
    if concretetype == 0 || concretetype == 1
        fc = param.mander.fc;
        Ec = param.mander.Ec;
        ecu = param.mander.ecu;
    elseif concretetype == 2
        SIGtc = param.hung.SIGtc;
        SIGtp = param.hung.SIGtp;
        EPSILONtc = param.hung.EPSILONtc;
        EPSILONtp = param.hung.EPSILONtp;
        EPSILONtu = param.hung.EPSILONtu;
        SIGcp = param.hung.SIGcp;
        SIGcu = param.hung.SIGcu;
        EPSILONcp = param.hung.EPSILONcp;
        EPSILONcu = param.hung.EPSILONcu;
    end
    
    ssAs = steelParam.ssAs;
    steelnum = steelParam.steelnum;
    totalAs = steelParam.totalAs;
    
    unr = visualParam.mander.unr;
    unecc = visualParam.mander.unecc;
    unfcc = visualParam.mander.unfcc;
    r = visualParam.mander.r;
    ecc = visualParam.mander.ecc;
    fcc = visualParam.mander.fcc;
    ectension = visualParam.mander.ectension;
    unconfined = visualParam.mander.unconfined;
    
    path = strcat(pathname, 'Mocur_', filename);
    fp = fopen(path, 'w');
    fprintf(fp, 'Mn                        Curvature                        c                    MaxStrainofCoreConc                  fsmax                 fsmin\r\n');
    
    % special case
    Paxial = P * b * h * fc;
    epsilonC = 0.00001;
    if Paxial > 0
        [epsilonC, xCurvature, yMn] = specialCase(param, steelParam, visualParam, Paxial, epsilonC);
    end

    PepsilonC = epsilonC;

    % write data to file
    c=h/2;central=h/2;flag=0;flagn=0;MaxStrainofCoreConc=0;Mpeak=0;
    v=2;epsilonTmax=0.001;Mnominal=0;Curvaturebefore=0;flagfcr=0;
    while (-MaxStrainofCoreConc<=ecu)&(epsilonTmax<=esu(sstyp(dtnum)))
        Csum=-1;Tsum=0.1;count=1;Mxx=0;Myy=0;
        c=central;
        while abs((-Csum/Tsum)-1)>0.01&&count<80
            Mcx=0;Mcy=0;Cc=0;Tc=0;Cs=0;T=0;Csum=0;Tsum=0;
            if concretetype==0||concretetype==1
                for g=1:1:40000%%%
                    epsilonc(g)=epsilonC*(crotation(2,g)-c)/c;   %斷面下壓上拉
                    if epsilonc(g)<0 && epsilonc(g)> -ecu    %受壓混凝土
                        if crotation(1,g)<=cover || crotation(1,g)>= (b-cover) || crotation(2,g)<=cover || crotation(2,g)>=(h-cover) %保護層
                            if epsilonc(g)>-0.004
                                eachfc(g)=-(unr*(-epsilonc(g)/unecc))/((unr-1)+(-epsilonc(g)/unecc).^unr)*unfcc;
                            elseif epsilonc(g)<=-0.004 && epsilonc(g)>=-0.006
                                eachfc(g)=-(-unconfined/0.002*-epsilonc(g)+3*unconfined);
                            else
                                eachfc(g)=0;
                            end
                        else
                            eachfc(g)=-(r*(-epsilonc(g)/ecc))/((r-1)+(-epsilonc(g)/ecc).^r)*fcc;  %將epsilonc(g)加負號變回正值，算完後再將應力變回負值代表壓力
                        end
                    elseif epsilonc(g)>0 && epsilonc(g)<=ectension       %受拉混凝土
                        eachfc(g)=epsilonc(g)*Ec;                        %受拉應力為正
                    else
                        eachfc(g)=0;
                    end
                    if  eachfc(g)<0
                        Cc=Cc+eachfc(g)*b/200*h/200;%%%
                    else
                        Tc=Tc+eachfc(g)*b/200*h/200;%%%
                    end

                    Mcx=-(eachfc(g)*b/200*h/200*((h/2)-cxy(2,g)))+Mcx;%%%
                    Mcy=(eachfc(g)*b/200*h/200*((b/2)-cxy(1,g)))+Mcy;%%%
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% HUNG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if concretetype==2
                for g=1:1:40000%%%
                    epsilonc(g)=epsilonC*(crotation(2,g)-c)/c;   %斷面下壓上拉
                    if epsilonc(g)>= 0 && epsilonc(g)<= EPSILONtc                 %受拉混凝土
                        if  EPSILONtc==0||SIGtc==0
                            eachfc(g)=0;
                        else
                            eachfc(g)=epsilonc(g)/EPSILONtc*SIGtc;
                        end
                    elseif epsilonc(g)>=EPSILONtc && epsilonc(g)<=EPSILONtp       %受拉混凝土
                        if  EPSILONtp==0||EPSILONtc==0||SIGtc==0||SIGtp==0
                            eachfc(g)=0;
                        else
                            eachfc(g)=SIGtc+(SIGtp-SIGtc)*(epsilonc(g)-EPSILONtc)/(EPSILONtp-EPSILONtc);
                        end
                    elseif epsilonc(g)>=EPSILONtp && epsilonc(g)<= EPSILONtu       %受拉混凝土
                        if  EPSILONtu==0||EPSILONtp==0||SIGtp==0
                            eachfc(g)=0;
                        else
                            eachfc(g)=SIGtp*(1-(epsilonc(g)-EPSILONtp)/(EPSILONtu-EPSILONtp));
                        end

                    elseif epsilonc(g)>=EPSILONtu                                 %受拉混凝土
                        eachfc(g)=0;
                    elseif epsilonc(g)>=-EPSILONcp && epsilonc(g)< 0              %受壓混凝土
                        eachfc(g)=-SIGcp*(2*(epsilonc(g)/(-EPSILONcp))-(epsilonc(g)/(-EPSILONcp))^2);
                    elseif epsilonc(g)>=-EPSILONcu && epsilonc(g)< -EPSILONcp       %受壓混凝土
                        eachfc(g)=-SIGcp*(1-((epsilonc(g)-(-EPSILONcp))/((-EPSILONcu)-(-EPSILONcp)))*(1-SIGcu/SIGcp));
                    else
                        eachfc(g)=-SIGcu;
                    end
                    if  eachfc(g)<0
                        Cc=Cc+eachfc(g)*b/200*h/200;%%%
                    else
                        Tc=Tc+eachfc(g)*b/200*h/200;%%%
                    end

                    Mcx=-(eachfc(g)*b/200*h/200*((h/2)-cxy(2,g)))+Mcx;%%%
                    Mcy=(eachfc(g)*b/200*h/200*((b/2)-cxy(1,g)))+Mcy;%%%
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%   Steel   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            epsilonTmax=epsilonC*(srotation(2,dtnum)-c)/c;


            Msx=0;Msy=0;Mn=0;fsmax=-20000;fsmin=20000;
            for m=1:1:steelnum
                epsilon(m)=epsilonC*(srotation(2,m)-c)/c;
                ety=fy(sstyp(m))/Es(sstyp(m));
                if epsilon(m)<ety&& epsilon(m)>=0
                    fs(m)=epsilon(m)*Es(sstyp(m));
                elseif epsilon(m)<0&& epsilon(m)>-(ety)&&flagfcr==0
                    fs(m)=epsilon(m)*Es(sstyp(m));
                elseif epsilon(m)<=esh(sstyp(m))&&epsilon(m)>=ety
                    fs(m)=fy(sstyp(m));
                elseif epsilon(m)>=-esh(sstyp(m))&&epsilon(m)<=-(ety)
                    if fcr(sstyp(m))<=fy(sstyp(m))%&&epsilon(m)>=0.004%%%%%%%%%%%%%%%%%%%%%%%%%%
                        fs(m)=0;
                        flagfcr=flagfcr+1;
                    elseif flagfcr==0
                        fs(m)=-fy(sstyp(m));
                    end
                elseif epsilon(m)<-esh(sstyp(m))&&epsilon(m)>=-esu(sstyp(m))
                    fs(m)=-(fu(sstyp(m))+(fy(sstyp(m))-fu(sstyp(m)))*((esu(sstyp(m))+epsilon(m))/(esu(sstyp(m))-esh(sstyp(m)))).^power(sstyp(m)));
                    if fs(m)<=-fcr(sstyp(m))
                        fs(m)=0;
                        flagfcr=flagfcr+1;
                    end
                elseif epsilon(m)>esh(sstyp(m))&&epsilon(m)<=esu(sstyp(m))
                    fs(m)=fu(sstyp(m))+(fy(sstyp(m))-fu(sstyp(m)))*((esu(sstyp(m))-epsilon(m))/(esu(sstyp(m))-esh(sstyp(m)))).^power(sstyp(m));
                else
                    fs(m)=0;
                end

                if fs(m)>fsmax
                    fsmax=fs(m);
                end
                if fs(m)<fsmin
                    fsmin=fs(m);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                epsilonc(m)=epsilonC*(srotation(2,m)-c)/c;
                if concretetype==0||concretetype==1
                    if epsilonc(m)<0 && epsilonc(m)> -ecu
                        eachfc(m)=-(r*(-epsilonc(m)/ecc))/((r-1)+(-epsilonc(m)/ecc).^r)*fcc;
                    end
                elseif concretetype==2
                    if epsilonc(m)>=-EPSILONcp && epsilonc(m)< 0              %受壓混凝土
                        eachfc(m)=-SIGcp*(2*(epsilonc(m)/(-EPSILONcp))-(epsilonc(m)/(-EPSILONcp))^2);
                    elseif epsilonc(m)>=-EPSILONcu && epsilonc(m)< -EPSILONcp       %受壓混凝土
                        eachfc(m)=-SIGcp*(1-((epsilonc(m)-(-EPSILONcp))/((-EPSILONcu)-(-EPSILONcp)))*(1-SIGcu/SIGcp));
                    else
                        eachfc(m)=-SIGcu;
                    end
                else
                    eachfc(m)=0;
                end
                fs(m)=fs(m)-eachfc(m);        %eachfc(m)為負值 fs(m)為負值


                if  fs(m)<0
                    Cs=Cs+fs(m)*ssAs(1,m);
                else
                    T=T+fs(m)*ssAs(1,m);
                end
                Msx=-fs(m)*ssAs(1,m)*((h/2)-srotation(2,m))+Msx;
                Msy=fs(m)*ssAs(1,m)*((b/2)-srotation(1,m))+Msy;
            end
            Csum=Cc+Cs;
            Tsum=T+Tc+Paxial;       %%%%%%%%%%%%%%%%%Paxial
            Mxx=(Mcx+Msx);
            Myy=(Mcy+Msy);

            %%%%%%%%%fprintf(fp,'%2.4f           %e         %f\r\n',c,epsilonC,-Csum/Tsum);%%%%%%%%%%%
            if c>3*h
                if (abs((-Csum/Tsum)-1))>0.01
                    c=c/(-Csum/Tsum);
                end
            else
                if (abs((-Csum/Tsum)-1))>0.01
                    c=c/(((-Csum/Tsum)+5)/6);
                end
            end
            count=count+1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        central=c;
    %     MaxStrainofCoreConc=-epsilonC/c*(c-syfirst(1,1));
        MaxStrainofCoreConc=-epsilonC/c*(c-(h - dtmax));
        Curvature=epsilonC/c;
        MM=(Mxx*Mxx)+(Myy*Myy);
        Mn=sqrt(MM);
        epsilonTy=fy(sstyp(dtnum))/Es(sstyp(dtnum));


        if epsilonTmax<=esu(sstyp(dtnum)) & Curvature>=Curvaturebefore
            %%%%%%%%%fprintf(fp,'Mn                        Curvature                        c                    MaxStrainofCoreConc                  fsmax                 fsmin\r\n');
            fprintf(fp,'%e              %e                   %2.4f                     %+2.7f                       %+5.1f                %+5.1f\r\n',Mn,Curvature,c,MaxStrainofCoreConc,fsmax,fsmin);
            xCurvature(1,v)=Curvature;
            yMn(1,v)=Mn;
            if epsilonTmax>epsilonTy&&flag==0
                Myield=Mn;
                Curvatureyield=Curvature;
                cyield=c;
                flag=flag+1;
            end

            if epsilonC>=0.003&&flagn==0&&PepsilonC<0.003
                Mnominal=Mn;
                Curvaturenominal=Curvature;
                cnominal=c;
                flagn=flagn+1;
            end
            if Mn>Mpeak
                Mpeak=Mn;Curvaturepeak=Curvature;cpeak=c;
            end
            Mult=Mn;Curvatureult=Curvature;cult=c;
            v=v+1;Curvaturebefore=Curvature;
        end

        %%epsilonC=epsilonC+0.0005;
        if epsilonC<=0.003
            epsilonC=epsilonC+0.00003;
        elseif 0.003<epsilonC<=0.01
            epsilonC=epsilonC+0.0002;
        else
            epsilonC=epsilonC+0.0005;
        end

    end
    
    setHandles(handles, flag, flagn, ...
                Curvatureyield, Myield, cyield, Mnominal, Curvaturenominal, cnominal, ...
                Mpeak, Curvaturepeak, cpeak, Mult, Curvatureult, cult, ...
                xCurvature, yMn);
end

function [epsilonC, xCurvature, yMn] = specialCase(param, steelParam, visualParam, Paxial, epsilonC)
    % load parameters
    concretetype = param.concretetype;
    sstyp = param.sstyp;
        
    h = param.h;
    b = param.b;
    fy = param.fy;
    Es = param.Es;
    
    if concretetype == 2
        SIGcp = param.hung.SIGcp;
        EPSILONcp = param.hung.EPSILONcp;
    end
    
    ssAs = steelParam.ssAs;
    steelnum = steelParam.steelnum;
    totalAs = steelParam.totalAs;
    
    unr = visualParam.mander.unr;
    unecc = visualParam.mander.unecc;
    unfcc = visualParam.mander.unfcc;
    
    Csum = -1;
%     Cs = 0;
    count = 1;
%     fsmax = -20000;
%     fsmin = 20000;
    Curvature = 0;
    
    while abs((-Csum / Paxial) - 1) > 0.03 && count < 200
%         Csum = 0;
        Cs = 0;
        Mpo = 0;
        fs(1:steelnum) = 0;
        for m = 1:steelnum
            if epsilonC * Es(sstyp(m)) >= fy(sstyp(m))
                fs(m) = -(fy(sstyp(m)));
                Cs = Cs + ssAs(1,m) * fs(m);
            else
                fs(m) = (-epsilonC) * Es(sstyp(m));
                Cs = Cs + ssAs(1,m) * fs(m);
            end
            Mpo = -(ssAs(1,m) * fs(m) * ((h / 2) - srotation(2, m))) + Mpo;
        end
        
        if concretetype == 0 || concretetype == 1
            eachconcreteforce = -(unr * (epsilonC / unecc)) / ((unr - 1) + (epsilonC / unecc) .^ unr) * unfcc;
        elseif concretetype == 2
            eachconcreteforce = -SIGcp * (2 * (-epsilonC / (-EPSILONcp)) - (-epsilonC / (-EPSILONcp)) ^ 2);
        end

        Csum = eachconcreteforce * b * h + Cs + totalAs * -eachconcreteforce;
        if abs((-Csum / Paxial) - 1) > 0.3 && (-Csum / Paxial) < 1
            epsilonC = epsilonC + 0.0001;
        elseif abs((-Csum / Paxial) - 1) > 0.1 && (-Csum / Paxial) < 1
            epsilonC = epsilonC + 0.00005;
        elseif abs((-Csum / Paxial) - 1) > 0.01 && (-Csum / Paxial) < 1
            epsilonC = epsilonC + 0.00002;
        elseif abs((-Csum / Paxial) - 1) > 0.01 && (-Csum / Paxial) > 1
            epsilonC = epsilonC - 0.00001;
        end
        count = count + 1;
    end
    fsmax = fs(1);
    fsmin = fs(1);
%     PepsilonC=epsilonC;
    c = Inf;
    MaxStrainofCoreConc = -epsilonC;
    xCurvature(1,1) = Curvature;
    yMn(1,1) = Mpo;
    fprintf(fp, '%e              %e                   %2.4f                          %+2.7f                       %+5.1f                %+5.1f\r\n', ...
        Mpo, Curvature, c, MaxStrainofCoreConc, fsmax, fsmin);
end

function setHandles(handles, flag, flagn, ...
                Curvatureyield, Myield, cyield, Mnominal, Curvaturenominal, cnominal, ...
                Mpeak, Curvaturepeak, cpeak, Mult, Curvatureult, cult, ...
                xCurvature, yMn)
    if flag == 0
        set(handles.text16, 'String', 'NA');
        set(handles.text17, 'String', 'NA');
        set(handles.text18, 'String', 'NA');
    else
        set(handles.text16, 'String', num2str(Curvatureyield, '%.3e'));
        set(handles.text17, 'String', num2str(Myield, '%.3e'));
        set(handles.text18, 'String', num2str(cyield, '%.2f'));
    end

    if flagn == 0
        set(handles.text19, 'String', 'NA');
        set(handles.text20, 'String', 'NA');
        set(handles.text21, 'String', 'NA');
    else
        set(handles.text19, 'String', num2str(Mnominal, '%.3e'));
        set(handles.text20, 'String', num2str(Curvaturenominal, '%.3e'));
        set(handles.text21, 'String', num2str(cnominal, '%.2f'));
    end

    set(handles.text31, 'String', num2str(Mpeak, '%.3e'));
    set(handles.text30, 'String', num2str(Curvaturepeak, '%.3e'));
    set(handles.text32, 'String', num2str(cpeak, '%.2f'));

    set(handles.text23, 'String', num2str(Mult, '%.3e'));
    set(handles.text24, 'String', num2str(Curvatureult, '%.3e'));
    set(handles.text25, 'String', num2str(cult, '%.2f'));

    axes(handles.axes4)
    plot(xCurvature, yMn, '-o', 'MarkerSize', 5);
    legend_list = handles.legend_list;
    legend(legend_list);
end