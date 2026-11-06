r"""core / hodge_grr."""

from sage.all import *

from functools import lru_cache
from math import factorial
from log_glsm.core.hodge import (
    HodgeIntegralBackend,
    HodgeIntegralRequest,
)

__all__ = ['GenusFourHodgeIntegrals', 'mumford_monomial']


@lru_cache(maxsize=None)
def mumford_monomial(genus, indices):
    """Expand a Hodge monomial in square-free lambda monomials.

    c(E)c(E^vee)=1 gives lambda_i^2 =
    2 sum_{k=1}^i (-1)^(k+1) lambda_{i-k} lambda_{i+k}.
    Every replacement strictly increases the sum of squared indices;
    hence this recursion terminates. Lambda_0=1, lambda_j=0 for j>g.
    """
    indices=tuple(sorted(i for i in indices if i))
    if any(i<0 or i>genus for i in indices):return ()
    for i in indices:
        if indices.count(i)<2:continue
        rest=list(indices);rest.remove(i);rest.remove(i);answer={}
        for k in range(1,i+1):
            if i+k>genus:continue
            term=tuple(sorted(rest+([i-k] if i-k else [])+[i+k]))
            for key,co in mumford_monomial(genus,term):
                answer[key]=answer.get(key,QQ.zero())+2*(-1)**(k+1)*co
        return tuple((key,co) for key,co in sorted(answer.items()) if co)
    # Hodge classes pull back from Mbar_g (Mbar_1,1 in genus one).
    if sum(indices)>(3*genus-3 if genus>=2 else genus):return ()
    return ((indices,QQ.one()),)

class GenusFourHodgeIntegrals(HodgeIntegralBackend):
    """Numerical Mumford GRR recursion; no boundary-class expansion.

    ch_k(E) = B_(k+1)/(k+1)! (kappa_k - sum psi_i^k
      + 1/2 xi_* sum_{a=0}^{k-1} (-1)^a psi_+^a psi_-^(k-1-a)).
    Odd Chern characters are additive under separating clutching and pull
    back unchanged under nonseparating clutching. See the source audit.
    """
    __hash__=object.__hash__
    @lru_cache(maxsize=None)
    def character_integral(self,genus,powers,characters):
        n=len(powers)
        if 2*genus-2+n<=0 or sum(powers)+sum(characters)!=3*genus-3+n:return QQ.zero()
        if not characters:return self.psi.integral(genus,powers)
        if genus==0:return QQ.zero()
        k=characters[-1];rest=characters[:-1]
        I=lambda g,p,c:self.character_integral(g,tuple(sorted(p,reverse=True)),tuple(sorted(c)))
        value=I(genus,powers+(k+1,),rest)
        for i in range(n):
            ps=list(powers);ps[i]+=k;value-=I(genus,tuple(ps),rest)
        boundary=QQ.zero()
        # Partitions carry their multiplicity, including repeated marks/chars.
        partitions={}
        for mask in range(1<<n):
            left=tuple(powers[i] for i in range(n) if mask&(1<<i))
            right=tuple(powers[i] for i in range(n) if not mask&(1<<i))
            partitions[left,right]=partitions.get((left,right),0)+1
        chpartitions={}
        for mask in range(1<<len(rest)):
            left=tuple(rest[i] for i in range(len(rest)) if mask&(1<<i))
            right=tuple(rest[i] for i in range(len(rest)) if not mask&(1<<i))
            chpartitions[left,right]=chpartitions.get((left,right),0)+1
        for a in range(k):
            term=I(genus-1,powers+(a,k-1-a),rest)
            for (lp,rp),pm in partitions.items():
                for (lc,rc),cm in chpartitions.items():
                    for h in range(genus+1):
                        # Dimension detects the only possible genus split.
                        if sum(lp)+a+sum(lc)!=3*h-2+len(lp):continue
                        term+=pm*cm*I(h,lp+(a,),lc)*I(genus-h,rp+(k-1-a,),rc)
            boundary+=(-1)**a*term
        return QQ(bernoulli(k+1))/factorial(k+1)*(value+boundary/2)

    @lru_cache(maxsize=None)
    def _lambda_characters(self,indices):
        P=PolynomialRing(QQ,2,('ch1','ch3'));a,b=P.gens()
        lam=(P.one(),a,a*a/2,a**3/6+2*b,a**4/24+2*a*b)
        polynomial=prod(lam[i] for i in indices)
        return tuple((tuple([1]*powers[0]+[3]*powers[1]),co) for powers,co in polynomial.dict().items())

    def integral(self,request_or_genus,psi_powers=None,lambda_indices=()):
        request=(request_or_genus if isinstance(request_or_genus,HodgeIntegralRequest)
                 else HodgeIntegralRequest(request_or_genus,psi_powers or (),lambda_indices))
        cached=self._lookup(request)
        if cached is not None:return cached
        if request.genus>4:raise ValueError('the Chern-character conversion is bounded at genus four')
        n=len(request.psi_powers);g=request.genus
        if sum(request.psi_powers)+sum(request.lambda_indices)!=3*g-3+n or 2*g-2+n<=0:return self._remember(request,QQ.zero())
        reduced=mumford_monomial(g,request.lambda_indices)
        answer=QQ.zero()
        for indices,scale in reduced:
            if not indices or indices==(g,) or indices==self.top_hodge_triple_indices(g):
                answer+=scale*super().integral(g,request.psi_powers,indices)
            else:
                for chars,co in self._lambda_characters(indices):
                    answer+=scale*co*self.character_integral(g,request.psi_powers,chars)
        return self._remember(request,answer)
