classdef KetySchmidtModel
    
    properties
        
        art_fit = 0;
        ven_fit = 0;
        lambda  = 1;
    end
    
    methods
        
        function this = KetySchmidtModel(afit, vfit)
            
            this.art_fit = afit;
            this.ven_fit = vfit;
        end
        
        function a1 = a(this, vasc)
            if (nargin < 2); vasc = 'v'; end
            if (strncmp('arter', vasc, 1))
                a1 = coeffvalues(this.art_fit);
            else
                a1 = coeffvalues(this.ven_fit);
            end
            a1 = a1(1);
        end
        
        function c1 = c(this, vasc)
            if (nargin < 2); vasc = 'v'; end
            if (strncmp('arter', vasc, 1))
                c1 = coeffvalues(this.art_fit);
            else
                c1 = coeffvalues(this.ven_fit);
            end
            c1 = c1(3);
        end
        
        function b1 = b(this, vasc)
            if (nargin < 2); vasc = 'v'; end
            if (strncmp('arter', vasc, 1))
                b1 = coeffvalues(this.art_fit);
            else
                b1 = coeffvalues(this.ven_fit);
            end
            b1 = b1(2);
        end
        
        function d1 = d(this, vasc)
            if (nargin < 2); vasc = 'v'; end
            if (strncmp('arter', vasc, 1))
                d1 = coeffvalues(this.art_fit);
            else
                d1 = coeffvalues(this.ven_fit);
            end
            d1 = d1(4);
        end
        
        function t = t0(this, vasc)
            if (nargin < 2); vasc = 'v'; end
            t = this.d(vasc) - log(this.c(vasc)/this.a(vasc) + 1)/this.b(vasc);
        end
        
        function t = tinf(this, tol)
            if (nargin < 2); tol = 0.0001; end
            t = 1;
            while (feval(this.ven_fit, inf) - feval(this.ven_fit, t) > tol)
                t = t*10;
            end
        end
        
        function F1 = F(this)
            
            num  = feval(this.ven_fit, inf)/this.lambda
            %t1   = this.t0('a'):1:this.tinf;
            %t2   = this.t0('v'):1:this.tinf;
            inta = integrate(this.art_fit, this.tinf, this.t0('a'));
            intv = integrate(this.ven_fit, this.tinf, this.t0('v'));
            F1   = num/(inta - intv);
        end
    end
end
