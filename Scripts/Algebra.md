Os comandos a seguir devem ser utilizados no software wxMaxima para o desenvolvimento algébrico da técnica de Freudenstein.

```
kill(all);

eq1: l2*cos(alpha)=l0-l1*cos(theta)+l3*cos(phi); 

eq2: l2*sin(alpha)=-l1*sin(theta)+l3*sin(phi);

eq3: eq1^2+eq2^2; 

expand(eq3);

eqt:trigsimp(eq3);

eqt:expand(%);

dir: rhs(eqt)$
esq: lhs(eqt)$
eqr: dir-esq;

k1t:coeff(eqr, sin(phi));
k2t:coeff(eqr, cos(phi));
k3t: ratsimp(eqr-((k1t*sin(phi))+k2t*cos(phi)));

eqteste: l0^2+l1^2+l3^2-l2^2-2*l0*l1*cos(theta)+2*l3*(l0-l1*cos(theta))*cos(phi)-2*l1*l3*sin(theta)*sin(phi);
ratsimp(eqr-eqteste);

seno: 2*t/(1+t^2)$
coseno: (1-t^2)/(1+t^2)$

eq4: k1*sin(phi)+k2*cos(phi)+k3;
eq4: subst(seno, sin(phi), eq4)$
eq4: subst(coseno, cos(phi), eq4)$
ratsimp(eq4);

r1: solve([eq4], [t]);
r2: r1[1];
r3t: subst(tan(phi/2),t,r2);
r3: subst(tan(phi/2),t,r1);
```
