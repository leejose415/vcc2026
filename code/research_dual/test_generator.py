"""Constraint tests, including moment conservation after narrowing and integer rounding."""
import numpy as np
from generator import DualGenerator,jackknife_by_gene

def test_constraints():
    rng=np.random.default_rng(61)
    template=rng.gamma(20,1,(400,40));template/=template.sum(1,keepdims=True)
    p=template.mean(0);q=p.copy();var=template.var(0)
    called=np.arange(40)%2==0
    engine=DualGenerator(p,q,var,template)
    results=[]
    for m in [1.,.5,.25]:
        X,stats=engine.generate(np.ones(40),called,m)
        assert np.issubdtype(X.dtype,np.integer) and X.min()>=0
        assert np.all(X.sum(1)<=1000000)
        assert stats['percell_max_relative_error']<.01
        assert stats['pooled_max_relative_error']<.01
        Y=X/X.sum(1)[:,None]
        assert np.all(Y.var(0)<=var*(1+1e-8)+1e-20)
        results.append((X,Y))
    for X,Y in results[1:]:
        assert np.array_equal(X.sum(0),results[0][0].sum(0))
        np.testing.assert_allclose(Y.mean(0),results[0][1].mean(0),atol=1e-14)
    assert results[2][1][:,called].var(0).sum()<results[1][1][:,called].var(0).sum()<results[0][1][:,called].var(0).sum()
    X=results[0][0];S=X.sum(0);L=X.sum(1)
    brute=np.log1p(50000*(S[None]-X)/(L.sum()-L[:,None]))
    expected=(len(X)-1)/len(X)*((brute-brute.mean(0))**2).sum(0)
    np.testing.assert_allclose(jackknife_by_gene(X),expected,rtol=1e-10,atol=1e-14)
    print('Passed: both integer moments, variance ceiling, narrowing, and exact jackknife contributions.')

if __name__=='__main__':test_constraints()
