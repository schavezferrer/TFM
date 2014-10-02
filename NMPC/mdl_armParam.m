

clear L
deg = pi/180;

L(1) = Revolute('d', -0.02625, 'a', 0.0301, 'alpha',pi/2, ...
    'I', [0, 0.35, 0, 0, 0, 0], ... 
    'r', [0, 0, 0], ... //ok
    'm', 0.02, ... 
    'Jm', 0, ... //ok
    'G', 0, ...
    'B', 0, ... //ok
    'Tc', [0 0], ...
    'qlim', [-67.5 67.5]*deg);

L(2) = Revolute('d', 0, 'a', 0.079647103, 'alpha', 0, ...
    'I', [0.13, 0.524, 0.539, 0, 0, 0], ...
    'r', [1, 0, 0], ...
    'm', 0.06, ... 
    'Jm', 0, ... // ok
    'G', 0, ...
    'B', 0, ...
    'Tc', [0 0], ...
    'qlim', [0 135]*deg, ... //ok
    'offset', pi+0.325); 

L(3) = Revolute('d', 0, 'a', 0.01581, 'alpha', pi/2,  ...
    'I', [0.066, 0.086, 0.0125, 0, 0, 0], ...
    'r', [-0.0203, -0.0141, 0.070], ...
    'm', 0.02, ...
    'Jm', 200e-6, ...
    'G', -53.7063, ...
    'B', 0, ...
    'Tc', [0 0], ...
    'qlim', [0 115]*deg, ... //ok
    'offset', pi/2-0.325);

L(4) = Revolute('d', -0.118, 'a', 0, 'alpha', pi,  ...
    'I', [1.8e-3, 1.3e-3, 1.8e-3, 0, 0, 0], ...
    'r', [0, 0.019, 0], ...
    'm', 0.02, ...
    'Jm', 33e-6, ...
    'G', 0, ...
    'B', 0, ...
    'Tc', [0 0], ...
    'qlim', [0 170]*deg);

L(5) = Revolute('d', 0, 'a', 0, 'alpha', -pi/2,  ...
    'I', [0.3e-3, 0.4e-3, 0.3e-3, 0, 0, 0], ...
    'r', [0, 0, 0], ...
    'm', 0.02, ...
    'Jm', 33e-6, ...
    'G', 0, ...
    'B', 0, ...
    'Tc', [0 0], ...
    'qlim', [-30 30]*deg);

L(6) = Revolute('d', 0, 'a', 0, 'alpha', 0,  ...
    'I', [0.3e-3, 0.4e-3, 0.3e-3, 0, 0, 0], ...
    'r', [0, 0, 0], ...
    'm', 0.02, ...
    'Jm', 33e-6, ...
    'G', 0, ...
    'B', 0, ...
    'Tc', [0 0], ...
    'qlim', [-85 85]*deg,...
    'offset', pi/2);


arm = SerialLink(L, 'name', 'Arm', ...
    'manufacturer', 'IRI', 'comment', 'Left Parameters');

arm.payload(0,[0 0 0]);

qn = [0 0 0 0 0 0];

clear L








