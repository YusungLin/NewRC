function param = refactor_read_data(handles)
    % Read s.
    %   read s from file and check for the s validation.
    
    [filename, pathname] = uigetfile({'*.txt'},'File Selector');
    fullpathname = strcat(pathname,filename);
    inputFile = fopen(fullpathname,'r');
    
    % 讀取開頭參數
    % read the parameters from the header lines
    [name, contype] = textread(fullpathname,'%s %s',1,'emptyvalue', NaN);
    handles.legend_list = [handles.legend_list; string(name)];
    
    concretetype = getConcreteType(contype);
    if concretetype == 0 || concretetype == 1
        [fc, Ec, K, ecu, ft] = textread(fullpathname,'%f %f %f %f %f',1,'headerlines',1,'emptyvalue', NaN);
    elseif concretetype == 2
        [SIGtc, SIGtp, EPSILONtc, EPSILONtp, EPSILONtu, SIGcp, SIGcu, EPSILONcp, EPSILONcu] = textread(fullpathname,'%f %f %f %f %f %f %f %f %f',1,'headerlines',1);
        fc = SIGcp;
        ecu = EPSILONcu;
    end
    
    [h, b, cover, P, stnum] = textread(fullpathname,'%f %f %f %f %d', 1, 'headerlines', 2);
    [fy, fu, Es, esh, esu, power, fcr] = textread(fullpathname, '%f %f %f %f %f %f %f', stnum, 'headerlines', 3, 'emptyvalue', NaN);
    
    % 檢查資料合法性
    % check for the s validation
    [fc, Ec, K, ecu, P, ...
        fy, fu, Es, esh, power, fcr] = headerValidation(concretetype, fc, Ec, K, ecu, P, ...
                                                           fy, fu, Es, esh, power, fcr);

    % discard the first stnum + 3 lines
    for num = 1:(stnum + 3)
        fgetl(inputFile);
    end
    % 讀取後續檔案
    % read the s from the rest of the file
    n = 1;
    while ~feof(inputFile)
        % Q: break or continue
%         l = fgetl(fr);
%         w=str2num(l);
%         if ~isnumeric(w)
%             break;
%         end
        dataRead = fgetl(inputFile);
        
%         count=0;kt=1;
%         while kt<=length(w)&& w(kt)>0
%             s(n,kt)=w(kt);
%             count=count+1;
%             kt=kt+1;
%         end
%         number(n)=count;
%         n=n+1;
        s(n, :) = str2double(regexp(dataRead, ' *', 'split'));
        tempSize = size(s(n,:));
        number(n) = tempSize(2);
        n = n + 1;
    end
    mat = size(s);
    
    dataValidation(h, s);

%     for i=1:1:mat(1)
%         for j=3:1:mat(2)
%             ss(i,j-2)=s(i,j);
%         end
%     end
    ss = s(:, 3:end);

%     g=1;
%     for n1=1:1:mat(1)
%         for n2=1:1:mat(2)-2
%             if ss(n1,n2)>0
%                 for num=1:1:stnum
%                     if s(n1,2)==fy(num)
%                         sstyp(1,g)=num;
%                         g=g+1;
%                     end
%                 end
%             elseif ss(n1,n2)==0
%                 sstyp(1,g)=0;
%                 g=g+1;
%             end
%         end
%     end
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     p=1;
%     for i=1:1:mat(1)
%         for j=1:1:mat(2)-2
%             ssrow(1,p)=[ss(i,j)];
%             p=p+1;
%         end
%     end
%     ssrow(find(ssrow==0))=[];
%     sstyp(find(sstyp==0))=[];
    for i = 1:mat(1)
        sstyp(i, 1:(mat(2) - 2)) = find(fy == s(i, 2));
    end
    sstyp(ss == 0) = 0;
    sstyp = reshape(sstyp', 1, []);
    ssrow = reshape(ss', 1, []);
    
    sstyp(sstyp == 0) = [];
    ssrow(ssrow == 0) = [];
    
    % 打包參數
    % pack all the parameters
    param.handles = handles;
    
    param.pathname = pathname;
    param.filename = filename;
    param.name = name;
    param.concretetype = concretetype;
        
    param.h = h;
    param.b = b;
    param.cover = cover;
    param.P = P;
    param.stnum = stnum;
    
    param.fy = fy;
    param.fu = fu;
    param.Es = Es;
    param.esh = esh;
    param.esu = esu;
    param.power = power;
    param.fcr = fcr;

    param.s = s;
    param.ss = ss;
    param.number = number;
    param.mat = mat;
    param.sstyp = sstyp;
    param.ssrow = ssrow;
    
    if concretetype == 0 || concretetype == 1
        mander.fc = fc;
        mander.Ec = Ec;
        mander.K = K;
        mander.ecu = ecu;
        mander.ft = ft;
        param.mander = mander;
    elseif concretetype == 2
        hung.SIGtc = SIGtc;
        hung.SIGtp = SIGtp;
        hung.EPSILONtc = EPSILONtc;
        hung.EPSILONtp = EPSILONtp;
        hung.EPSILONtu = EPSILONtu;
        hung.SIGcp = SIGcp;
        hung.SIGcu = SIGcu;
        hung.EPSILONcp = EPSILONcp;
        hung.EPSILONcu = EPSILONcu;
        param.hung = hung;
    end
end

function concretetype = getConcreteType(contype)
    if strcmp(contype,'NaN') == 1
        concretetype = 0;
    elseif strcmp(contype,'MANDER') == 1
        concretetype = 1;
    elseif strcmp(contype,'HUNG') == 1
        concretetype = 2;
    end
end

function promptError(msg)
    handler = msgbox(msg,'輸入錯誤','error');
    set(handler,'color','w');
    waitfor(handler);
    close(gcf);
end

function [fc, Ec, K, ecu, P, ...
    fy, fu, Es, esh, power, fcr] ...
    = headerValidation(concretetype, fc, Ec, K, ecu, P, ...
    fy, fu, Es, esh, power, fcr)

    if concretetype == 0 || concretetype == 1
        if fc > 1500 || fc < 100
            ('fc 輸入錯誤');
        end
        if isnan(Ec)
            Ec = 15000*sqrt(fc);
        end
        if K < 1
            promptError('K 輸入錯誤');
        end
        if isnan(K)
            K = 1;
        end
        if isnan(ecu)
            ecu = 0.004;
        end
    end

%     for num=1:1:stnum
%         if isnan(fcr(num))
%             fcr(num)=fu(num)+100;
%         end
%     end
    fcr(isnan(fcr)) = fu(isnan(fcr)) + 100;
    
    if P > 0.99
        promptError('P 輸入錯誤');
    end
    
%     for num=1:1:stnum
%         if fy(num)>10000||fy(num)<1000
%             hm=msgbox(['fy輸入錯誤'],'輸入錯誤','error');
%             set(hm,'color','w');
%             waitfor(hm);
%             close(gcf);
%         end
%         if fu(num)<fy(num)
%             hm=msgbox(['fu輸入錯誤'],'輸入錯誤','error');
%             set(hm,'color','w');
%             waitfor(hm);
%             close(gcf);
%         end
%         if esh(num)<(fy(num)/Es(num))
%             hm=msgbox(['esh輸入錯誤'],'輸入錯誤','error');
%             set(hm,'color','w');
%             waitfor(hm);
%             close(gcf);
%         end
%         if power(num)<1
%             hm=msgbox(['power輸入錯誤'],'輸入錯誤','error');
%             set(hm,'color','w');
%             waitfor(hm);
%             close(gcf);
%         end
%     end
    if sum(fy > 10000 | fy < 1000) > 0
        promptError('fy 輸入錯誤');
    end
    if sum(fu < fy) > 0
        promptError('fu 輸入錯誤');
    end
    if sum(esh < fy ./ Es) > 0
        promptError('esh 輸入錯誤');
    end
    if sum(power < 1) > 0
        promptError('power 輸入錯誤');
    end
end

function dataValidation(h, s)
%     for n=1:1:mat(1)
%         if s(n,1)>h
%             hm=msgbox(['di輸入錯誤'],'輸入錯誤','error');
%             set(hm,'color','w');
%             waitfor(hm);
%             close(gcf);
%         end
%         if s(n,2)>10000||s(n,2)<1000
%             hm=msgbox(['fy輸入錯誤'],'輸入錯誤','error');
%             set(hm,'color','w');
%             waitfor(hm);
%             close(gcf);
%         end
% 
%         for q=3:1:mat(2)
%             re=rem(s(n,q),1);
%             if s(n,q)>=19||s(n,q)<0||s(n,q)==1||s(n,q)==2||s(n,q)==13||s(n,q)==15||s(n,q)==17||re~=0   %re~=0  代表輸入非整數
%                 hm=msgbox(['鋼筋號數輸入錯誤'],'輸入錯誤','error');
%                 set(hm,'color','w');
%                 waitfor(hm);
%                 close(gcf);
%             end
%         end
%     end
    if sum(s(:, 1) > h) > 0
        promptError('di 輸入錯誤');
    end
    if sum(s(:, 2) > 10000 | s(:, 2) < 1000) > 0
        promptError('fy 輸入錯誤');
    end
    if sum(sum(s(:, 3:end) >= 19 | s(:, 3:end) <= 0 | ...
            s(:, 3:end) == 1 | s(:, 3:end) == 2 | s(:, 3:end) == 13 | s(:, 3:end) == 15 | s(:, 3:end) == 17 | ...
            rem(s(:, 3:end), 1) ~= 0)) > 0  % re ~= 0 代表輸入非整數 means non-integer
        promptError('鋼筋號數輸入錯誤');
    end
end