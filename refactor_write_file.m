function refactor_write_file(param, steelParam, visualParam, ...
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
    c = h / 2;
%     central = h / 2;
    flag = 0;
    flagn = 0;
    MaxStrainofCoreConc = 0;
    Mpeak = 0;
    v = 2;
    epsilonTmax = 0.001;
    Mnominal = 0;
    Curvaturebefore = 0;
    while (-MaxStrainofCoreConc <= ecu) && (epsilonTmax <= esu(sstyp(dtnum)))
%         c = central;
        
        [Mxx, Myy, c, fsmax, fsmin] = refactor_data_preprocess(param, steelParam, visualParam, ...
                                epsilonC, c, cxy, srotation, crotation, Paxial);
        epsilonTmax = epsilonC * (srotation(2, dtnum) - c) / c;
        
    %     MaxStrainofCoreConc=-epsilonC/c*(c-syfirst(1,1));
        MaxStrainofCoreConc = -epsilonC / c * (c - (h - dtmax));
        Curvature = epsilonC / c;
%         MM=(Mxx*Mxx)+(Myy*Myy);
%         Mn=sqrt(MM);
        Mn = sqrt(Mxx * Mxx + Myy * Myy);
        epsilonTy = fy(sstyp(dtnum)) / Es(sstyp(dtnum));


        if epsilonTmax <= esu(sstyp(dtnum)) && Curvature >= Curvaturebefore
            %%%%%%%%%fprintf(fp,'Mn                        Curvature                        c                    MaxStrainofCoreConc                  fsmax                 fsmin\r\n');
            fprintf(fp,'%e              %e                   %2.4f                     %+2.7f                       %+5.1f                %+5.1f\r\n', ...
                Mn, Curvature, c, MaxStrainofCoreConc, fsmax, fsmin);
            xCurvature(1,v) = Curvature;
            yMn(1,v) = Mn;
            if epsilonTmax > epsilonTy && flag == 0
                Myield = Mn;
                Curvatureyield = Curvature;
                cyield = c;
                flag = flag + 1;
            end

            if epsilonC >= 0.003 && flagn == 0 && PepsilonC < 0.003
                Mnominal = Mn;
                Curvaturenominal = Curvature;
                cnominal = c;
                flagn = flagn + 1;
            end
            if Mn > Mpeak
                Mpeak = Mn;
                Curvaturepeak = Curvature;
                cpeak = c;
            end
            Mult = Mn;
            Curvatureult = Curvature;
            cult = c;
            v = v + 1;
            Curvaturebefore = Curvature;
        end

        %%epsilonC=epsilonC+0.0005;
        if epsilonC <= 0.003
            epsilonC = epsilonC + 0.00003;
        elseif 0.003 < epsilonC && epsilonC <= 0.01
            epsilonC = epsilonC + 0.0002;
        else
            epsilonC = epsilonC + 0.0005;
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
        % preallocation
        fs = zeros(1, steelnum);
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