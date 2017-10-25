classdef TestDWT3 < matlab.unittest.TestCase
    %
    % USAGE:
    % >> testCase = TestDWT3;
    % >> res = run(testCase)
    %
    properties
        X, N, eps, nres;
    end  
    
    methods (Test)
        % Constructor
        function obj = TestDWT3(testCase)
            obj.N = 2^6;
            obj.X = rand([obj.N, obj.N, obj.N]);
            obj.eps = 1e-5;
            obj.nres = 2;
        end
        
        function obj = testDecRecQuadratic(testCase)

            N = testCase.N;
            X = testCase.X;
            nres = testCase.nres;
            eps = testCase.eps;

            Z = DWT3Impl(X, nres, 'pwl0', 'per');
            Y = IDWT3Impl(Z, nres, 'pwl0', 'per');

            err = norm1(Y-X);
            testCase.verifyTrue(err < eps);

        end

        function obj = testDecRecRectangular(testCase)

            N = testCase.N;
            M = N*2;
            X = rand(M,N,N);
            nres = testCase.nres;
            eps = testCase.eps;

            Z = DWT3Impl(X, nres, 'cdf97', 'sym');
            Y = IDWT3Impl(Z, nres, 'cdf97', 'sym');

            err = norm1(Y-X);
            testCase.verifyTrue(err < eps);

        end


        function test_BD_precond_vm2(testCase)
            testPrecondBoundary(testCase, 2);
        end
        function test_BD_precond_vm3(testCase)
            testPrecondBoundary(testCase, 3);
        end
        function test_BD_precond_vm4(testCase)
            testPrecondBoundary(testCase, 4);
        end
        function test_BD_precond_vm5(testCase)
            testPrecondBoundary(testCase, 5);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       Daubechies wavelets                            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testDB2(testCase)
            testWaveletDecRec(testCase, 'db2');
        end
        function testDB3(testCase)
            testWaveletDecRec(testCase, 'db3');
        end
        function testDB4(testCase)
            testWaveletDecRec(testCase, 'db4');
        end
        
        function testDB2bd(testCase)
            testWaveletDecRecBoundary(testCase, 'db2');
        end
        function testDB3bd(testCase)
            testWaveletDecRecBoundary(testCase, 'db3');
        end
        function testDB4bd(testCase)
            testWaveletDecRecBoundary(testCase, 'db4');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       Symlets wavelets                            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testSYM2(testCase)
            testWaveletDecRec(testCase, 'sym2');
        end
        function testSYM3(testCase)
            testWaveletDecRec(testCase, 'sym3');
        end
        function testSYM4(testCase)
            testWaveletDecRec(testCase, 'sym4');
        end
        
        function testSYM2bd(testCase)
            testWaveletDecRecBoundary(testCase, 'sym2');
        end
        function testSYM3bd(testCase)
            testWaveletDecRecBoundary(testCase, 'sym3');
        end
        function testSYM4bd(testCase)
            testWaveletDecRecBoundary(testCase, 'sym4');
        end

    end
    methods (Access=private)

        
        function testPrecondBoundary(testCase, vm)
            
            wave_name = sprintf('db%d', vm);
            
            eps = 1e-5; 
            nres = 2;
            N = 64;

            
            idx = linspace(0,1,N);
            [X, Y, Z] = meshgrid(idx, idx, idx);
            A = X+Y+Z;
            
            idxl = 1:N/2;
            idxh = N/2+1:N;
            
            success = 1;
            
            for i = 1:vm
                R = A.^(i-1);
                T = DWT3Impl(R, nres, wave_name, 'bd_pre');
                a1 = norm1(T(idxl,idxh,idxh));
                a2 = norm1(T(idxl,idxl,idxh));
                a3 = norm1(T(idxh,idxl,idxh));
                a4 = norm1(T(idxh,idxh,idxh));
                a5 = norm1(T(idxl,idxh,idxl));
                a6 = norm1(T(idxh,idxl,idxl));
                a7 = norm1(T(idxh,idxh,idxl));
                %a1+a2+a3+a4+a5+a6+a7
                success = success & (a1+a2+a3+a4+a5+a6+a7) < eps;
            end
            testCase.verifyTrue(success);
            
        end

        function testWaveletDecRec(testCase, wave_name)
            X = testCase.X;
            eps = testCase.eps;
            nres = testCase.nres;
            Z = IDWT3Impl(DWT3Impl(X, nres, wave_name, 'per', 0), nres, wave_name, 'per', 0);
            testCase.verifyTrue(norm1(Z-X) < eps);
            Z = IDWT3Impl(DWT3Impl(X, nres, wave_name, 'symm', 0), nres, wave_name, 'symm', 0);
            testCase.verifyTrue(norm1(Z-X) < eps);
            Z = IDWT3Impl(DWT3Impl(X, nres, wave_name, 'per', 1), nres, wave_name, 'per', 1);
            testCase.verifyTrue(norm1(Z-X) < eps);
            Z = IDWT3Impl(DWT3Impl(X, nres, wave_name, 'symm', 1), nres, wave_name, 'symm', 1);
            testCase.verifyTrue(norm1(Z-X) < eps);
        end
        
        function testWaveletDecRecBoundary(testCase, wname)
            X = testCase.X;
            eps = testCase.eps;
            nres = testCase.nres;
            
            Z = IDWT3Impl(DWT3Impl(X, nres, wname, 'bd', 0), nres, wname, 'bd', 0);
            testCase.verifyTrue(norm1(Z-X) < eps);
            
            Z = IDWT3Impl(DWT3Impl(X, nres, wname, 'bd_pre', 0), nres, wname, 'bd_pre', 0);
            testCase.verifyTrue(norm1(Z-X) < eps);
        end
        
    end

end 


function y = norm1(X)
    y = sum(abs(X(:)));
end

