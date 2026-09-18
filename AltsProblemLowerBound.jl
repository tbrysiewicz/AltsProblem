using HomotopyContinuation

function AltsProblem()
    @var x a y b xhat ahat yhat bhat
    @var γ[1:8] γhat[1:8] δ[1:8] δhat[1:8]

    D1 = [(ahat*x - δhat[i]*x)*γ[i] + (a*xhat - δ[i]*xhat)*γhat[i] +
          (ahat - xhat)*δ[i] + (a - x)*δhat[i] - δ[i]*δhat[i] for i in 1:8]
    D2 = [(bhat*y - δhat[i]*y)*γ[i] + (b*yhat - δ[i]*yhat)*γhat[i] +
          (bhat - yhat)*δ[i] + (b - y)*δhat[i] - δ[i]*δhat[i] for i in 1:8]
    D3 = [γ[i]*γhat[i] + γ[i] + γhat[i] for i in 1:8]

    System(vcat(D1, D2, D3);
           variables=[x, a, y, b, xhat, ahat, yhat, bhat, γ..., γhat...],
           parameters=[δ..., δhat...])
end
