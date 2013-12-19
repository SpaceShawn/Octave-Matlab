% requires a differential equation containing two depenedent variables
function  [x_next y_next]=SB_RungeKutta(Fun,x_i,y_i,h)
    K1 = Fun(x_i,y_i);
    K2 = Fun(x_i+h, y_i+h);
    x_next = x_i + h;
    y_next = y_i + (1/2)*(K1+K2)*h;
end;

x = 0:0.1:2;
input_function = -1.2*y+7*e.^(-0.3*x);

next_values  = SB_RungeKutta(Fun,0,3,0.5)
