classdef KetySchmidt
    
    properties
        
        Vdead   = 1.2/1000; % The dead-space for the arterial catheter and valve was 0.2 mL
                            % for the venous catheter:  1.0 mL
                            % for blood-collection syringes:  0.03 mL
        Vdsyr   = 0.0; % Rinsed syringes with heplock
        Vtube   = 10.0/1000;
        Vdetect = 333.0/1000; % tubing, gas analyzer, etc.
        Vsyr    = [];
        ppm     = [];
        epochs  = [];
        c       = [];
        Nunits    = 'mol';
        Vunits    = 'L';
        concUnits = 'M';
        Avogadro  = 6.02214179e23; % #/mole
    end
        
    methods
        
        function this = KetySchmidt(ep, vs, ppm)
            
            %% KETYSCHMIDT ctor
            %  cf. Lee et al., AJNR 2010
            %  Usage: this = KetySchmidt(ep, vs, ppm)
            %                            ^   ^   ^ col-vec
            %                                ^ in mL
            %
            this.epochs = ep;
            this.Vsyr   = vs - this.Vdsyr;
            this.Vsyr   = this.Vsyr/1000;
            this.ppm    = ppm;
            assert(length(ep) == length(vs));
            assert(length(vs) == length(ppm));
        end
        
        function g1 = g(this, vs) % vector of samples
            
            %% G is the volume fraction of blood left behind 
            %  after syringe extraction at time-point m.
            %  Usage:  g1 = obj.g(vs) 
            %                     ^
            %        vs is a vector
            %
            % $\forall m | g_m = frac{mod(V_{s,m}, V_d)}{V_d}$
            %
            if (nargin < 2); vs = this.Vsyr; end
            g1 = mod(vs, this.Vdead) ./ this.Vdead;
        end
        
        function c1 = get.c(this)
            
            %% GET.C retrieves this.c with lazy initialization
            %  units of M sample, dead-space-corrected for in vivo
            %
            if ~isempty(this.c)
                c1 = this.c;
            else
                c1    = zeros(length(this.Vsyr), 1);
                c_iv  = this.c_invivo;
                c1(1) = c_iv(1);                % M
                Ns    = this.Nsyringe;          % mol
                assert(length(this.Vsyr) > 2);
                for m = 2:length(this.Vsyr)
                    gm    = this.g(m)*this.g(m-1);
                    c1(m) = Ns(m) + gm*this.Vdead*c1(m-1); % mol
                    c1(m) = c1(m) / (this.Vsyr(m) + gm*this.Vdead); 
                end
            end
            c1 = this.conc2torr(c1);
            assert(isnumeric(c1));
        end
        
        function c1 = c_invivo(this)
            
            %% C_INVIVO in M
            %
            c1 = this.Nsyringe;
            c1 = c1 ./ this.Vsyr;
        end
        
        function n1 = Nsyringe(this)
            
            %% Nsyringe 
            %  units specified by this.Nunits
            %
            n1 = (this.ppm / 1e6) * this.Nchamber;
        end
        
        function n1 = Nchamber(this)
            
            %% NCHAMBER n = P*V/R*T
            %  units specified by this.Nunits
            %
            V  = this.Vdetect + this.Vtube;
            assert(strcmp(this.Vunits, 'L'));
            T  = 298;      % K
            P  = 760;      % Torr
            R  = 62.36367; % L Torr K^-1 mol^-1
            n1 = P*V/(R*T);
        end
        
        function c1 = conc2torr(this, c0)
            
            %% CONC2TORR P = (n/V)*R*T
            %
            T  = 298;      % K
            R  = 62.36367; % Torr K^-1 M^-1
            assert(strcmp('M', this.concUnits));
            c1 = c0*R*T;
        end
        
        function n1 = convert(this, lbl, n0)
            
            %% CONVERT expects moles, rescales units per this.Nunits
            %
            switch (lbl)
                case 'moles'
                    switch (this.Nunits)
                        case {'mol','moles'}
                            n1 = n0;
                        case 'mmol'
                            n1 = 1.0e3*n0;
                        case 'micromol'
                            n1 = 1.0e6*n0;
                        case 'nmol'
                            n1 = 1.0e9*n0;
                        case 'pmol'
                            n1 = 1.0e12*n0;
                        case 'fmol'
                            n1 = 1.0e15*n0;
                        otherwise
                            error('mlperf:NotImplemented', ['KetySchmidt.Nunits->' this.Nunits]);
                    end
                case 'L'
                    switch (this.Vunits)
                        case 'L'
                            n1 = n0;
                        case 'mL'
                            n1 = 1.0e3*n0;
                        case 'microL'
                            n1 = 1.0e6*n0;
                        case 'nanoL'
                            n1 = 1.0e9*n0;
                        case 'pL'
                            n1 = 1.0e12*n0;
                        case 'fL'
                            n1 = 1.0e15*n0;
                        otherwise
                            error('mlperf:NotImplemented', ['KetySchmidt.Vunits->' this.Vunits]);
                    end
                case 'molar'
                    switch (this.concUnits)
                        case {'molar','M'}
                            n1 = n0;
                        case {'mmolar','mM'}
                            n1 = 1.0e3*n0;
                        case {'micromolar','microM'}
                            n1 = 1.0e6*n0;
                        case {'nmolar','nM'}
                            n1 = 1.0e9*n0;
                        case {'pmolar','pM'}
                            n1 = 1.0e12*n0;
                        case {'fmolar','fM'}
                            n1 = 1.0e15*n0;
                        otherwise
                            error('mlperf:NotImplemented', ['KetySchmidt.concUnits->' this.concUnits]);
                    end
                otherwise
                    error('mlperf:NotImplemented', ['KetySchmidt.convert.lbl->' lbl]);
            end
            
        end
        
        function createfigure(X1, Y1, X2, Y2, X3, YMatrix1, Y3, Y4)
            %CREATEFIGURE(X1,Y1,X2,Y2,X3,YMATRIX1,Y3,Y4)
            %  X1:  vector of x data
            %  Y1:  vector of y data
            %  X2:  vector of x data
            %  Y2:  vector of y data
            %  X3:  vector of x data
            %  YMATRIX1:  matrix of y data
            %  Y3:  vector of y data
            %  Y4:  vector of y data
            
            %  Auto-generated by MATLAB on 20-Jan-2010 17:25:51
            
            % Create figure
            figure1 = figure('XVisual','0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)','Renderer','painters',...
                'Color',[1 1 1]);
            
            % Create axes
            axes1 = axes('Parent',figure1,...
                'Position',[0.127324749642346 0.49728555917481 0.777675250357654 0.42771444082519],...
                'FontWeight','bold',...
                'FontSize',14);
            % Uncomment the following line to preserve the X-limits of the axes
            % xlim(axes1,[5 30]);
            % Uncomment the following line to preserve the Y-limits of the axes
            % ylim(axes1,[0 80]);
            box(axes1,'on');
            hold(axes1,'all');
            
            % Create plot
            plot(X1,Y1,'Parent',axes1,'Marker','o','LineWidth',2,'LineStyle','none','Color',[1 0 0],...
                'DisplayName','c1a vs. t1a');
            
            % Create plot
            plot(X2,Y2,'Parent',axes1,'Marker','square','LineWidth',2,'LineStyle','none',...
                'Color',[0.47843137254902 0.0627450980392157 0.894117647058824],...
                'DisplayName','c1v vs. t1v');
            
            % Create multiple lines using matrix input to plot
            plot1 = plot(X3,YMatrix1,'Parent',axes1);
            set(plot1(1),'LineWidth',2,'Color',[0.266666666666667 0.0431372549019608 0.494117647058824],...
                'DisplayName','  fit 1v');
            set(plot1(2),'LineStyle',':','Color',[0.266666666666667 0.0431372549019608 0.494117647058824],...
                'DisplayName','    Pred bnds (fit 1v)');
            set(plot1(3),'LineWidth',2,'Color',[0.847058823529412 0.16078431372549 0],'DisplayName','  fit 1a');
            set(plot1(4),'LineStyle',':','Color',[0.847058823529412 0.16078431372549 0],...
                'DisplayName','    Pred bnds (fit 1a)');
            
            % Create xlabel
            xlabel('time / min','FontWeight','bold','FontSize',14);
            
            % Create ylabel
            ylabel('intravascular pNO_2 / mm Hg','FontWeight','bold','FontSize',14);
            
            % Create title
            title('NO_2 Tracer Kinetics 32 mm Hg ETCO_2','FontWeight','bold','FontSize',18);
            
            % Create axes
            axes2 = axes('Parent',figure1,'Position',[0.127324749642346 0.11 0.777675250357654 0.274364820846906],...
                'FontWeight','bold',...
                'FontSize',14);
            % Uncomment the following line to preserve the X-limits of the axes
            % xlim(axes2,[5 30]);
            % Uncomment the following line to preserve the Y-limits of the axes
            % ylim(axes2,[-15 15]);
            box(axes2,'on');
            hold(axes2,'all');
            
            % Create plot
            plot(X2,Y3,'Parent',axes2,'Marker','square','LineWidth',2,'LineStyle','none',...
                'Color',[0.266666666666667 0.0431372549019608 0.494117647058824],...
                'DisplayName','fit 1v');
            
            % Create plot
            plot(X1,Y4,'Parent',axes2,'Marker','o','LineWidth',2,'LineStyle','none',...
                'Color',[0.847058823529412 0.16078431372549 0],...
                'DisplayName','fit 1a');
            
            % Create title
            title('Residuals','FontWeight','bold','FontSize',18);
            
            % Create legend
            legend1 = legend(axes1,'show');
            set(legend1,'Interpreter','none','Visible','off',...
                'Position',[0.619888128471818 0.530817944920632 0.245350500715308 0.11957111834962],...
                'FontSize',11,...
                'FontWeight','normal');
            
            % Create legend
            legend2 = legend(axes2,'show');
            set(legend2,'Interpreter','none','Visible','off',...
                'Position',[0.719957081545065 0.144543973941368 0.136623748211731 0.0538816503800217],...
                'FontSize',11,...
                'FontWeight','normal');
            
        end
    end
    
    methods (Static)
        
        function ks = factory(filename)
            
            %% FACTORY
            import mlperf.*;
            dstruct = KetySchmidt.importfile(filename);
            ks.runs = cell(size(dstruct.runs));
            for r = 1:length(ks.runs)
                len   = length(dstruct.runs{r}.arterial);
                ep_   = zeros(1,len);
                Vsyr_ = zeros(1,len);
                ppm_  = zeros(1,len);
                for t = 1:len
                    ep_(t)   = dstruct.runs{r}.arterial{t}.minutes;
                    Vsyr_(t) = dstruct.runs{r}.arterial{t}.Vsyr;
                    ppm_(t)  = dstruct.runs{r}.arterial{t}.ppm;
                end
                ks.runs{r}.arterial = KetySchmidt(ep_, Vsyr_, ppm_);
                len   = length(dstruct.runs{r}.venous);
                ep_   = zeros(1,len);
                Vsyr_ = zeros(1,len);
                ppm_  = zeros(1,len);
                for t = 1:len
                    ep_(t)   = dstruct.runs{r}.venous{t}.minutes;
                    Vsyr_(t) = dstruct.runs{r}.venous{t}.Vsyr;
                    ppm_(t)  = dstruct.runs{r}.venous{t}.ppm;
                end
                ks.runs{r}.venous   = KetySchmidt(ep_, Vsyr_, ppm_);
            end
        end % factory
        
        function dstruct = importfile(filename)
            
            %% IMPORTFILE(FILENAME)
            %  Imports data from tab-delimited text files such as Lois/2006mar1/2006mar1.txt
            %  Expects:  
            %    sampleId	Vsyr	hr	min	sec	ppm
            %    arterial1
            %    a01	0.28	0	6	41	12
            %    ...
            %    venous1
            %    v01	0.28	0	7	15	1
            %    ...
            %
            %  FILETOREAD1:  tab-delimited text file to read

            %  Auto-generated by MATLAB on 19-Jan-2010 19:10:26

            % Import the file
            dstruct = importdata(filename);

            % Create new variables in the base workspace from those fields.
            vars = fieldnames(dstruct);
            for i = 1:length(vars)
                assignin('base', vars{i}, dstruct.(vars{i}));
            end
            
            % Check type & completenes of variables
            patterns = { 'sample', 'V', 'hr', 'min', 'sec', 'ppm' };
            for p = 1:length(patterns)
                assert(all(strfind(dstruct.textdata{1,p}, patterns{p})));
            end
            assert(all((size(dstruct.data) + [1 1]) == (size(dstruct.textdata))));
            
            % Determining data array sizes
            art_lbl_idx = strmatch('arterial', dstruct.textdata);
            art_dat_idx = art_lbl_idx - 1;
            ven_lbl_idx = strmatch('venous',   dstruct.textdata);
            ven_dat_idx = ven_lbl_idx - 1;
            len_labels  = size(dstruct.textdata,1);
            len_data    = size(dstruct.data,1);
            if (len_labels && (isempty(art_lbl_idx) || isempty(ven_lbl_idx)))
                error('mlperf:MissingInputErr', 'KetySchmidt.importfile needs arterial and venous label-lines');
            end
            assert(length(art_lbl_idx) == length(ven_lbl_idx));
            last_art   = ven_dat_idx - 1;
            N_art_sets = last_art - art_dat_idx;
            last_ven   = circshift(art_dat_idx,1) - 1; last_ven(end) = len_data;
            N_ven_sets = last_ven - ven_dat_idx;
            
            % Populating dstruct
            assert(length(N_art_sets) == length(N_ven_sets)); % every art run must have a ven run
            len_runs     = length(N_art_sets);
            dstruct.runs = cell(len_runs,1);
            for r = 1:len_runs

                sId   = cell(N_art_sets(r),1);
                for s = 1:N_art_sets(r)
                    sId{s,1} = dstruct.textdata{art_lbl_idx(r)+s,1};
                end
                Vsyr_ = dstruct.data(    art_dat_idx(r)+1:last_art(r),1);
                hr    = dstruct.data(    art_dat_idx(r)+1:last_art(r),2);
                mn    = dstruct.data(    art_dat_idx(r)+1:last_art(r),3);
                sc    = dstruct.data(    art_dat_idx(r)+1:last_art(r),4);
                mn    = mn + 60*hr + sc/60;
                ppm_  = dstruct.data(    art_dat_idx(r)+1:last_art(r),5);
                for s = 1:N_art_sets(r)
                    dstruct.runs{r}.arterial{s} = struct('sampleId', sId{s}, 'Vsyr', Vsyr_(s), 'minutes', mn(s), 'ppm', ppm_(s));
                end
                
                sId   = cell(N_ven_sets(r),1);
                for s = 1:N_ven_sets(r)
                    sId{s,1} = dstruct.textdata{ven_lbl_idx(r)+s,1};
                end
                Vsyr_ = dstruct.data(    ven_dat_idx(r)+1:last_ven(r),1);
                hr    = dstruct.data(    ven_dat_idx(r)+1:last_ven(r),2);
                mn    = dstruct.data(    ven_dat_idx(r)+1:last_ven(r),3);
                sc    = dstruct.data(    ven_dat_idx(r)+1:last_ven(r),4);
                mn    = mn + 60*hr + sc/60;
                ppm_  = dstruct.data(    ven_dat_idx(r)+1:last_ven(r),5);
                for s = 1:N_ven_sets(r)
                    dstruct.runs{r}.venous{s}   = struct('sampleId', sId{s}, 'Vsyr', Vsyr_(s), 'minutes', mn(s), 'ppm', ppm_(s));
                end
            end
        end % importfile
    end
end
