classdef TestDecRec < matlab.unittest.TestCase
    % TestDecRec test if a the DWT and IDWT function calles reveres each other
    % for all kernel functions. 
    %
    % USAGE:
    % >> testCase = TestDecRec;
    % >> res = run(testCase)
    %
    properties
        x, N, eps, nres;
    end  
    
    methods (Test)
        % Constructor
        function obj = TestDecRec(testCase)
            obj.N = 2^7;
            obj.x = rand([obj.N,1]);
            obj.eps = 1e-8;
            obj.nres = 2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function testCDF97(testCase)
            testWaveletDecRec(testCase, 'cdf97');
        end
        function testCDF53(testCase)
            testWaveletDecRec(testCase, 'spline53');
        end
        function testpwl0(testCase)
            testWaveletDecRec(testCase, 'pwl0');
        end
        function testpwl2(testCase)
            testWaveletDecRec(testCase, 'pwl2');
        end
        function testHaar(testCase)
            testWaveletDecRec(testCase, 'haar');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       Daubechies wavelets                            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This function fails
        %function testDB1(testCase)
        %    testWaveletDecRec(testCase, 'db1');
        %end
        function testDB2(testCase)
            testWaveletDecRec(testCase, 'db2');
        end
        function testDB3(testCase)
            testWaveletDecRec(testCase, 'db3');
        end
        function testDB4(testCase)
            testWaveletDecRec(testCase, 'db4');
        end
        function testDB5(testCase)
            testWaveletDecRec(testCase, 'db5');
        end
        function testDB6(testCase)
            testWaveletDecRec(testCase, 'db6');
        end
        function testDB7(testCase)
            testWaveletDecRec(testCase, 'db7');
        end
        function testDB8(testCase)
            testWaveletDecRec(testCase, 'db8');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function testDB2bd(testCase)
            testWaveletDecRecBoundary(testCase, 'db2');
        end
        function testDB3bd(testCase)
            testWaveletDecRecBoundary(testCase, 'db3');
        end
        function testDB4bd(testCase)
            testWaveletDecRecBoundary(testCase, 'db4');
        end
        function testDB5bd(testCase)
            testWaveletDecRecBoundary(testCase, 'db5');
        end
        function testDB6bd(testCase)
            testWaveletDecRecBoundary(testCase, 'db6');
        end
        function testDB7bd(testCase)
            testWaveletDecRecBoundary(testCase, 'db7');
        end
        %function testDB8bd(testCase)
        %    testWaveletDecRecBoundary(testCase, 'db8');
        %end
                
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
        function testSYM5(testCase)
            testWaveletDecRec(testCase, 'sym5');
        end
        function testSYM6(testCase)
            testWaveletDecRec(testCase, 'sym6');
        end
        function testSYM7(testCase)
            testWaveletDecRec(testCase, 'sym7');
        end
        function testSYM8(testCase)
            testWaveletDecRec(testCase, 'sym8');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function testSYM2bd(testCase)
            testWaveletDecRecBoundary(testCase, 'sym2');
        end
        function testSYM3bd(testCase)
            testWaveletDecRecBoundary(testCase, 'sym3');
        end
        function testSYM4bd(testCase)
            testWaveletDecRecBoundary(testCase, 'sym4');
        end
        function testSYM5bd(testCase)
            testWaveletDecRecBoundary(testCase, 'sym5');
        end
        %function testSYM6bd(testCase)
        %    testWaveletDecRecBoundary(testCase, 'sym6');
        %end
        %function testSYM7bd(testCase)
        %    testWaveletDecRecBoundary(testCase, 'sym7');
        %end
        %function testSYM8bd(testCase)
        %    testWaveletDecRecBoundary(testCase, 'sym8');
        %end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       Precond wavelets                            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function test_BD_precond_vm2(testCase)
            testPrecondBoundary_daubechies(testCase, 2);
        end
        function test_BD_precond_vm3(testCase)
            testPrecondBoundary_daubechies(testCase, 3);
        end
        function test_BD_precond_vm4(testCase)
            testPrecondBoundary_daubechies(testCase, 4);
        end
        function test_BD_precond_vm5(testCase)
            testPrecondBoundary_daubechies(testCase, 5);
        end
        function test_BD_precond_vm6(testCase)
            testPrecondBoundary_daubechies(testCase, 6);
        end

        function test_2d_input(testCase)
            eps = testCase.eps;
            N = 32;
            nres = 2;
            wname = 'pwl0';
            bd_mode = 'per';

            X = rand(N,N);
            Z = dwt_impl(X, nres, wname, bd_mode);
            Y = idwt_impl(Z, nres, wname, bd_mode);
            err = norm(X-Y,'fro');
            testCase.verifyTrue(err < eps);
        end

    end
    methods (Access=private)

        function testWaveletDecRec(testCase, wave_name)
            x = testCase.x;
            eps = testCase.eps;
            nres = testCase.nres;

            z = idwt_impl(dwt_impl(x, wave_name, nres, 'per'), wave_name, nres, 'per');
            testCase.verifyTrue(norm(z-x,2) < eps);
            z = idwt_impl(dwt_impl(x, wave_name, nres, 'symm'), wave_name, nres, 'symm');
            testCase.verifyTrue(norm(z-x,2) < eps);

        end

        function testWaveletDecRecBoundary(testCase, wave_name)
            x = testCase.x;
            eps = testCase.eps;
            nres = testCase.nres;

            z = idwt_impl(dwt_impl(x, wave_name, nres, 'bd'), wave_name, nres, 'bd');
            testCase.verifyTrue(norm(z-x,2) < eps);

            z = idwt_impl(dwt_impl(x, wave_name, nres, 'bd', 'bd_pre'), wave_name, nres, 'bd', 'bd_pre');
            testCase.verifyTrue(norm(z-x,2) < eps);
        end
        
        function testPrecondBoundary_daubechies(testCase, vm)
            
            wave_name = sprintf('db%d', vm);
            
            eps = testCase.eps;
            nres = testCase.nres;
            N = testCase.N;

            x = linspace(0,1,N)';
            success = 1;
            for i = 1:vm
                y = x.^(i-1);
                z = dwt_impl(x, wave_name, nres, 'bd', 'bd_pre');
                a = z(N/2^nres+1:N);

                success = success & norm(a,2) < eps;
            end
            testCase.verifyTrue(success);
            
        end
        
        function testPrecondBoundary_symlet(testCase, vm)
            
            wave_name = sprintf('sym%d', vm);
            
            eps = testCase.eps;
            nres = testCase.nres;
            N = testCase.N;

            x = linspace(0,1,N)';
            success = 1;
            for i = 1:vm
                y = x.^(i-1);
                z = dwt_impl(x, wave_name, nres, 'bd', 'bd_pre');
                a = z(N/2^nres+1:N);

                success = success & norm(a,2) < eps;
            end
            testCase.verifyTrue(success);

        end
    end
end 



