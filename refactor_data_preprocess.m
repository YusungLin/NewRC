function [Mxx, Myy, c, fsmax, fsmin] = refactor_data_preprocess(param, steelParam, visualParam, ...
                                            epsilonC, c, cxy, srotation, crotation, Paxial)
    concretetype = param.concretetype;
    sstyp = param.sstyp;
        
    h = param.h;
    b = param.b;
    cover = param.cover;
    
    fy = param.fy;
    fu = param.fu;
    Es = param.Es;
    esh = param.esh;
    esu = param.esu;
    power = param.power;
    fcr = param.fcr;
    
    flagfcr = 0;
    
    if concretetype == 0 || concretetype == 1
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
    
    unr = visualParam.mander.unr;
    unecc = visualParam.mander.unecc;
    unfcc = visualParam.mander.unfcc;
    r = visualParam.mander.r;
    ecc = visualParam.mander.ecc;
    fcc = visualParam.mander.fcc;
    ectension = visualParam.mander.ectension;
    unconfined = visualParam.mander.unconfined;
    
    Csum = -1;
    Tsum = 0.1;
    count = 1;
    Mxx = 0;
    Myy = 0;

    while abs((-Csum / Tsum) - 1) > 0.01 && count < 80
        Mcx = 0;
        Mcy = 0;
        Cc = 0;
        Tc = 0;
        Cs = 0;
        T = 0;
%         Csum = 0;
%         Tsum = 0;
        % preallocation
        epsilonc = zeros(1, 40000);
        eachfc = zeros(1, 40000);
        if concretetype == 0 || concretetype == 1
            for g=1:40000%%%
                epsilonc(g) = epsilonC * (crotation(2, g) - c) / c;   % 斷面下壓上拉
                if epsilonc(g) < 0 && epsilonc(g) > -ecu    % 受壓混凝土
                    if crotation(1, g) <= cover || crotation(1, g) >= (b - cover) || crotation(2, g) <= cover || crotation(2, g) >= (h - cover) % 保護層
                        if epsilonc(g) > -0.004
                            eachfc(g) = -(unr * (-epsilonc(g) / unecc)) / ((unr - 1) + (-epsilonc(g) / unecc) .^ unr) * unfcc;
                        elseif epsilonc(g) <= -0.004 && epsilonc(g) >= -0.006
                            eachfc(g) = -(-unconfined / 0.002 * -epsilonc(g) + 3 * unconfined);
                        else
                            eachfc(g) = 0;
                        end
                    else
                        eachfc(g) = -(r * (-epsilonc(g) / ecc)) / ((r - 1) + (-epsilonc(g) / ecc) .^ r) * fcc;  % 將epsilonc(g)加負號變回正值，算完後再將應力變回負值代表壓力
                    end
                elseif epsilonc(g) > 0 && epsilonc(g) <= ectension       % 受拉混凝土
                    eachfc(g) = epsilonc(g) * Ec;                        % 受拉應力為正
                else
                    eachfc(g) = 0;
                end
                if  eachfc(g) < 0
                    Cc = Cc + eachfc(g) * b / 200 * h / 200;
                else
                    Tc = Tc + eachfc(g) * b / 200 * h / 200;
                end

                Mcx = -(eachfc(g) * b / 200 * h / 200 * ((h / 2) - cxy(2, g))) + Mcx;
                Mcy = (eachfc(g) * b / 200 * h / 200 * ((b / 2) - cxy(1,g))) + Mcy;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% HUNG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if concretetype == 2
            for g = 1:40000
                epsilonc(g) = epsilonC * (crotation(2, g) - c) / c;     % 斷面下壓上拉
                if epsilonc(g) >= 0 && epsilonc(g) <= EPSILONtc         % 受拉混凝土
                    if  EPSILONtc == 0 || SIGtc == 0
                        eachfc(g) = 0;
                    else
                        eachfc(g) = epsilonc(g) / EPSILONtc * SIGtc;
                    end
                elseif epsilonc(g) >= EPSILONtc && epsilonc(g) <= EPSILONtp     % 受拉混凝土
                    if EPSILONtp == 0 || EPSILONtc == 0 || SIGtc == 0 || SIGtp == 0
                        eachfc(g) = 0;
                    else
                        eachfc(g) = SIGtc + (SIGtp - SIGtc) * (epsilonc(g) - EPSILONtc) / (EPSILONtp - EPSILONtc);
                    end
                elseif epsilonc(g) >= EPSILONtp && epsilonc(g) <= EPSILONtu     % 受拉混凝土
                    if EPSILONtu == 0 || EPSILONtp == 0 || SIGtp == 0
                        eachfc(g) = 0;
                    else
                        eachfc(g) = SIGtp * (1 - (epsilonc(g) - EPSILONtp) / (EPSILONtu - EPSILONtp));
                    end
                elseif epsilonc(g) >= EPSILONtu                                 % 受拉混凝土
                    eachfc(g) = 0;
                elseif epsilonc(g) >= -EPSILONcp && epsilonc(g) < 0             % 受壓混凝土
                    eachfc(g) = -SIGcp * (2 * (epsilonc(g) / -EPSILONcp) - (epsilonc(g) / -EPSILONcp) ^ 2);
                elseif epsilonc(g) >= -EPSILONcu && epsilonc(g) < -EPSILONcp    % 受壓混凝土
                    eachfc(g) = -SIGcp * (1 - ((epsilonc(g) - (-EPSILONcp)) / (-EPSILONcu - (-EPSILONcp))) * (1 - SIGcu / SIGcp));
                else
                    eachfc(g) = -SIGcu;
                end
                if  eachfc(g)<0
                    Cc = Cc + eachfc(g) * b / 200 * h / 200;
                else
                    Tc = Tc + eachfc(g) * b / 200 * h / 200;
                end

                Mcx = -(eachfc(g) * b / 200 * h / 200 * ((h / 2) - cxy(2, g))) + Mcx;
                Mcy = (eachfc(g) * b / 200 * h / 200 * ((b / 2) - cxy(1, g))) + Mcy;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%   Steel   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         epsilonTmax=epsilonC*(srotation(2,dtnum)-c)/c;


        Msx = 0;
        Msy = 0;
%         Mn = 0;
        fsmax = -20000;
        fsmin = 20000;
        % preallocation
        epsilon = zeros(1, steelnum);
        fs = zeros(1, steelnum);
        for m = 1:steelnum
            epsilon(m) = epsilonC * (srotation(2, m) - c) / c;
            ety = fy(sstyp(m)) / Es(sstyp(m));
            if epsilon(m) < ety && epsilon(m) >= 0
                fs(m) = epsilon(m) * Es(sstyp(m));
            elseif epsilon(m) < 0 && epsilon(m) > -(ety) && flagfcr == 0
                fs(m) = epsilon(m) * Es(sstyp(m));
            elseif epsilon(m) <= esh(sstyp(m)) && epsilon(m) >= ety
                fs(m) = fy(sstyp(m));
            elseif epsilon(m) >= -esh(sstyp(m)) && epsilon(m) <= -ety
                if fcr(sstyp(m)) <= fy(sstyp(m)) %&&epsilon(m)>=0.004%%%%%%%%%%%%%%%%%%%%%%%%%%
                    fs(m) = 0;
                    flagfcr = flagfcr + 1;
                elseif flagfcr == 0
                    fs(m) = -fy(sstyp(m));
                end
            elseif epsilon(m) < -esh(sstyp(m)) && epsilon(m) >= -esu(sstyp(m))
                fs(m) = -(fu(sstyp(m)) + (fy(sstyp(m)) - fu(sstyp(m))) * ((esu(sstyp(m)) + epsilon(m)) / (esu(sstyp(m)) - esh(sstyp(m)))) .^ power(sstyp(m)));
                if fs(m) <= -fcr(sstyp(m))
                    fs(m) = 0;
                    flagfcr = flagfcr + 1;
                end
            elseif epsilon(m) > esh(sstyp(m)) && epsilon(m) <= esu(sstyp(m))
                fs(m) = fu(sstyp(m)) + (fy(sstyp(m)) - fu(sstyp(m))) * ((esu(sstyp(m)) - epsilon(m)) / (esu(sstyp(m)) - esh(sstyp(m)))) .^ power(sstyp(m));
            else
                fs(m) = 0;
            end

            if fs(m) > fsmax
                fsmax = fs(m);
            end
            if fs(m) < fsmin
                fsmin = fs(m);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            epsilonc(m) = epsilonC * (srotation(2, m) - c) / c;
            if concretetype == 0 || concretetype == 1
                if epsilonc(m) < 0 && epsilonc(m) > -ecu
                    eachfc(m) = -(r * (-epsilonc(m) / ecc)) / ((r - 1) + (-epsilonc(m) / ecc) .^ r) * fcc;
                else
                    eachfc(m) = 0;
                end
            elseif concretetype == 2
                if epsilonc(m) >= -EPSILONcp && epsilonc(m) < 0                 % 受壓混凝土
                    eachfc(m) = -SIGcp * (2 * (epsilonc(m) / -EPSILONcp) - (epsilonc(m) / -EPSILONcp) ^ 2);
                elseif epsilonc(m) >= -EPSILONcu && epsilonc(m) < -EPSILONcp    % 受壓混凝土
                    eachfc(m) = -SIGcp * (1 - ((epsilonc(m) - (-EPSILONcp)) / (-EPSILONcu - (-EPSILONcp))) * (1 - SIGcu / SIGcp));
                else
                    eachfc(m) = -SIGcu;
                end
            else
                eachfc(m) = 0;
            end
            fs(m) = fs(m) - eachfc(m);        % eachfc(m) 為負值 fs(m) 為負值


            if  fs(m) < 0
                Cs = Cs + fs(m) * ssAs(1, m);
            else
                T = T + fs(m) * ssAs(1, m);
            end
            Msx = -fs(m) * ssAs(1,m) * ((h / 2) - srotation(2, m)) + Msx;
            Msy = fs(m) * ssAs(1,m) * ((b / 2) - srotation(1,m)) + Msy;
        end
        Csum = Cc + Cs;
        Tsum = T + Tc + Paxial;       %%%%%%%%%%%%%%%%%Paxial
        Mxx = Mcx + Msx;
        Myy = Mcy + Msy;

        %%%%%%%%%fprintf(fp,'%2.4f           %e         %f\r\n',c,epsilonC,-Csum/Tsum);%%%%%%%%%%%
        if c > 3 * h
            if abs((-Csum / Tsum) - 1) > 0.01
                c = c / (-Csum / Tsum);
            end
        else
            if abs((-Csum / Tsum) - 1) > 0.01
                c = c / (((-Csum / Tsum) + 5) / 6);
            end
        end
        count = count + 1;
    end
end