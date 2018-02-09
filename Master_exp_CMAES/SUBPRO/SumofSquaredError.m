function z = SumofSquaredError(x_experimental,x)
    z = sum((x_experimental-x).^2);
end