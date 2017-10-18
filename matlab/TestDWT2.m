classdef TestDWT2 < matlab.unittest.TestCase
    %
    % USAGE:
    % >> testCase = TestDWT2;
    % >> res = run(testCase)
    %
    properties
        X, N, eps, nres;
    end  
    
    methods (Test)
        % Constructor
        function obj = TestDWT2(testCase)
            obj.N = 2^7;
            obj.X = rand([obj.N,obj.N]);
            obj.eps = 1e-8;
            obj.nres = 2;
        end
        
        function obj = testDecRecQuadratic(testCase)

            N = testCase.N;
            X = testCase.X;
            nres = testCase.nres;
            eps = testCase.eps;

            Z = DWT2Impl(X, nres, 'pwl0', 'per');
            Y = IDWT2Impl(Z, nres, 'pwl0', 'per');

            err = norm(Y-X, 'fro');
            testCase.verifyTrue(err < eps);

        end

        function obj = testDecRecRectangular(testCase)

            N = testCase.N;
            M = N*2;
            X = rand(M,N);
            nres = testCase.nres;
            eps = testCase.eps;

            Z = DWT2Impl(X, nres, 'cdf97', 'per');
            Y = IDWT2Impl(Z, nres, 'cdf97', 'per');

            err = norm(Y-X, 'fro');
            testCase.verifyTrue(err < eps);

        end

        function obj = testDecRec3D(testCase)

            N = testCase.N;
            M = N*2;
            d = 5;
            X = rand(M,N, d);
            nres = testCase.nres;
            eps = testCase.eps;

            Z = DWT2Impl(X, nres, 'db2', 'per');
            Y = IDWT2Impl(Z, nres, 'db2', 'per');
            
            err = 0;
            for i = 1:d
                err = err + norm(Y(:,:,i)-X(:,:,i), 'fro');
            end
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
        function test_BD_precond_vm6(testCase)
            testPrecondBoundary(testCase, 6);
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
            
            eps = testCase.eps;
            nres = 2;
            N = 64;

            
            idx = linspace(0,1,N);
            [X, Y] = meshgrid(idx, idx);
            Z = X+Y;
            
            success = 1;
            for i = 1:vm
                R = Z.^(i-1);
                T = DWT2Impl(R, nres, wave_name, 'bd_pre');
                a1 = T(N/2^nres+1:N, :);
                a2 = T(1:N/2^nres, N/2^nres+1:N);
                success = success & (norm(a1,2) < eps) & (norm(a2, 2) < eps);
            end
            testCase.verifyTrue(success);
            
        end

        function testWaveletDecRec(testCase, wave_name)
            X = testCase.X;
            eps = testCase.eps;
            nres = testCase.nres;
            Z = IDWT2Impl(DWT2Impl(X, nres, wave_name, 'per', 0), nres, wave_name, 'per', 0);
            testCase.verifyTrue(norm(Z-X,'fro') < eps);
            Z = IDWT2Impl(DWT2Impl(X, nres, wave_name, 'symm', 0), nres, wave_name, 'symm', 0);
            testCase.verifyTrue(norm(Z-X,'fro') < eps);
            Z = IDWT2Impl(DWT2Impl(X, nres, wave_name, 'per', 1), nres, wave_name, 'per', 1);
            testCase.verifyTrue(norm(Z-X,'fro') < eps);
            Z = IDWT2Impl(DWT2Impl(X, nres, wave_name, 'symm', 1), nres, wave_name, 'symm', 1);
            testCase.verifyTrue(norm(Z-X,'fro') < eps);
        end
        
        function testWaveletDecRecBoundary(testCase, wname)
            X = testCase.X;
            eps = testCase.eps;
            nres = testCase.nres;
            
            Z = IDWT2Impl(DWT2Impl(X, nres, wname, 'bd', 0), nres, wname, 'bd', 0);
            testCase.verifyTrue(norm(Z-X, 'fro') < eps);
            
            Z = IDWT2Impl(DWT2Impl(X, nres, wname, 'bd_pre', 0), nres, wname, 'bd_pre', 0);
            testCase.verifyTrue(norm(Z-X,'fro') < eps);
        end


    end

end 




