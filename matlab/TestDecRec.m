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
            testWaveletDecRec(testCase, 'cdf53');
        end
        function testpwl0(testCase)
            testWaveletDecRec(testCase, 'pwl0');
        end
        function testpwl2(testCase)
            testWaveletDecRec(testCase, 'pwl2');
        end
        function testHaar(testCase)
            testWaveletDecRec(testCase, 'pwl2');
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
        function testDB8bd(testCase)
            testWaveletDecRecBoundary(testCase, 'db8');
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
        function testSYM6bd(testCase)
            testWaveletDecRecBoundary(testCase, 'sym6');
        end
        function testSYM7bd(testCase)
            testWaveletDecRecBoundary(testCase, 'sym7');
        end
        function testSYM8bd(testCase)
            testWaveletDecRecBoundary(testCase, 'sym8');
        end
    end
    methods (Access=private)
        function testWaveletDecRec(testCase, wave_name)
            x = testCase.x;
            eps = testCase.eps;
            nres = testCase.nres;
            z = IDWTImpl(DWTImpl(x, nres, wave_name, 0, 0), nres, wave_name, 0, 0);
            testCase.verifyTrue(norm(z-x,2) < eps);
            z = IDWTImpl(DWTImpl(x, nres, wave_name, 1, 0), nres, wave_name, 1, 0);
            testCase.verifyTrue(norm(z-x,2) < eps);
            z = IDWTImpl(DWTImpl(x, nres, wave_name, 0, 1), nres, wave_name, 0, 1);
            testCase.verifyTrue(norm(z-x,2) < eps);
            z = IDWTImpl(DWTImpl(x, nres, wave_name, 1, 1), nres, wave_name, 1, 1);
            testCase.verifyTrue(norm(z-x,2) < eps);
        end
        function testWaveletDecRecBoundary(testCase, wave_name)
            x = testCase.x;
            eps = testCase.eps;
            nres = testCase.nres;
            z = IDWTImpl(DWTImpl(x, nres, wave_name, 2, 0), nres, wave_name, 2, 0);
            testCase.verifyTrue(norm(z-x,2) < eps);
            z = IDWTImpl(DWTImpl(x, nres, wave_name, 3, 0), nres, wave_name, 3, 0);
            testCase.verifyTrue(norm(z-x,2) < eps);
        end
    end
end 



