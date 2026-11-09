// 运行：node source/verify.cjs —— 从 index.html 抽取 <script id="core"> 的纯数学内核并做独立数值检查
const fs=require('fs'),path=require('path');
const html=fs.readFileSync(path.join(__dirname,'..','index.html'),'utf8');
const F=new Function(html.match(/<script id="core">([\s\S]*?)<\/script>/)[1]+'\nreturn FourierCore;')();
let pass=0,fail=0;const ok=(name,cond,info='')=>{if(cond){pass++;console.log('✓',name,info)}else{fail++;console.log('✗',name,info)}};
const PI=Math.PI,P=F.presets,N=60;
const maxRel=(arr,ref)=>{let e=0;for(let n=1;n<=N;n++)e=Math.max(e,Math.abs(arr[n]-ref(n)));return e};
let c=F.coeffs(P.saw,'full',N);ok('锯齿 x：bₙ = 2(−1)ⁿ⁺¹/n',maxRel(c.b,n=>2*(-1)**(n+1)/n)<1e-5&&maxRel(c.a,()=>0)<1e-9);
c=F.coeffs(P.parab,'full',N);ok('x²：a₀ = 2π²/3，aₙ = 4(−1)ⁿ/n²',Math.abs(c.a[0]-2*PI*PI/3)<1e-6&&maxRel(c.a,n=>4*(-1)**n/n**2)<1e-5);
c=F.coeffs(P.square,'full',N);ok('方波：bₙ = 4/(nπ)（n 奇），0（n 偶）',maxRel(c.b,n=>n%2?4/(n*PI):0)<1e-8);
c=F.coeffs(P.abs,'full',N);ok('|x|：a₀ = π，aₙ = −4/(πn²)（n 奇）',Math.abs(c.a[0]-PI)<1e-6&&maxRel(c.a,n=>n%2?-4/(PI*n*n):0)<1e-6);
c=F.coeffs(P.parab,'sine',N);ok('x² 的正弦级数：bₙ = (2/π)[(−1)ⁿ⁺¹π²/n + 2((−1)ⁿ−1)/n³]',maxRel(c.b,n=>2/PI*((-1)**(n+1)*PI*PI/n+2*((-1)**n-1)/n**3))<1e-5);
// e^{cos x} 的系数 = 2 Iₙ(1)
const I=n=>{let s=0,t=1/[...Array(n).keys()].reduce((p,k)=>p*(k+1),1)/2**n;for(let k=0;k<30;k++){s+=t;t*=0.25/((k+1)*(k+1+n))}return s};
c=F.coeffs(P.smooth,'full',12);let e7=0;for(let n=0;n<=12;n++)e7=Math.max(e7,Math.abs(c.a[n]-2*I(n)));ok('e^{cos x}：aₙ = 2Iₙ(1)（指数衰减）',e7<1e-10,`max 误差 ${e7.toExponential(2)}`);
// Gibbs
c=F.coeffs(P.square,'full',400);let mx=0;for(let i=1;i<=4000;i++)mx=Math.max(mx,F.partial(c,400,i*0.05/4000));const g=(mx-1)/2;ok('方波 N=400：超调 / 跳跃 → 0.0895',Math.abs(g-0.08949)<0.001,`= ${g.toFixed(5)}`);
let fm=0;for(let i=1;i<=4000;i++)fm=Math.max(fm,F.partial(c,50,i*PI/4000,true));ok('Fejér 平均 σ₅₀ 不超出 f 的取值范围（无超调）',fm<=1+1e-9,`max σ = ${fm.toFixed(6)}`);
ok('跳跃点处 S_N(0) = 0 = ½[f(0−)+f(0+)]',Math.abs(F.partial(c,200,0))<1e-12);
c=F.coeffs(P.abs,'full',400);ok('|x| 在连续点 x=1：|S₄₀₀(1) − 1| < 10⁻³',Math.abs(F.partial(c,400,1)-1)<1e-3);
// Parseval
ok('(1/π)∫x² 的 L² 范数平方 = 2π²/3',Math.abs(F.normSq(P.saw,'full')-2*PI*PI/3)<1e-6);
c=F.coeffs(P.square,'full',400);let s=0;for(let n=1;n<=400;n++)s+=c.b[n]**2;let tail=0;for(let n=401;n<2e6;n+=2)tail+=16/(n*n*PI*PI);ok('方波 Parseval：Σbₙ²（含尾部）= (1/π)∫1 = 2',Math.abs(s+tail-2)<1e-5,`和 = ${(s+tail).toFixed(7)}`);
// 延拓
ok('正弦级数 = 奇延拓，余弦级数 = 偶延拓',Math.abs(F.ext(P.parab,'sine',-1)+1)<1e-15&&Math.abs(F.ext(P.parab,'cos',-1)-1)<1e-15&&Math.abs(F.ext(P.parab,'full',1+2*PI)-1)<1e-12);
console.log(`\n${pass} 项通过，${fail} 项失败`);process.exit(fail?1:0);
