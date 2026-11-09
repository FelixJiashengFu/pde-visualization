// 运行：node source/verify.cjs —— 从 index.html 抽取 <script id="core"> 的纯数学内核并做独立数值检查
const fs=require('fs'),path=require('path');
const html=fs.readFileSync(path.join(__dirname,'..','index.html'),'utf8');
const code=html.match(/<script id="core">([\s\S]*?)<\/script>/)[1];
const H=new Function(code+'\nreturn HeatCore;')();
let pass=0,fail=0;const ok=(name,cond,info='')=>{if(cond){pass++;console.log('✓',name,info)}else{fail++;console.log('✗',name,info)}};
const PI=Math.PI;
// 1 误差函数
const ref=[[0.5,0.5204998778],[1,0.8427007929],[2,0.9953222650],[-1.3,-0.9340079449]];
ok('erf 与参考值一致（|误差|<2e-7）',ref.every(([x,v])=>Math.abs(H.erf(x)-v)<2e-7));
// 2 两模态初值的正弦系数
const cm=H.coeffs(H.presets.modes,12,'dirichlet');ok('sin πx + ½ sin 6πx 的系数为 b₁=1、b₆=0.5、其余为 0',[...cm].every((v,n)=>Math.abs(v-(n===1?1:n===6?0.5:0))<1e-9));
// 3 方块初值系数与闭式 bₙ = 2/(nπ)[cos(0.3nπ) − cos(0.7nπ)]
const cb=H.coeffs(H.presets.box,30,'dirichlet');let e3=0;for(let n=1;n<=30;n++)e3=Math.max(e3,Math.abs(cb[n]-2/(n*PI)*(Math.cos(0.3*n*PI)-Math.cos(0.7*n*PI))));ok('方块初值系数与闭式一致',e3<1e-5,`max 误差 ${e3.toExponential(2)}`);
// 4 边界条件
const ct=H.coeffs(H.presets.tri,150,'dirichlet');ok('Dirichlet 级数在 x=0、1 处为 0',[0.001,0.02,0.1].every(t=>Math.abs(H.evalAt(ct,'dirichlet',1,t,0))<1e-12&&Math.abs(H.evalAt(ct,'dirichlet',1,t,1))<1e-12));
const cn=H.coeffs(H.presets.tri,150,'neumann');const dx=1e-5;ok('Neumann 级数在端点 uₓ≈0',[0.01,0.05].every(t=>Math.abs(H.evalAt(cn,'neumann',1,t,dx)-H.evalAt(cn,'neumann',1,t,0))/dx<1e-3&&Math.abs(H.evalAt(cn,'neumann',1,t,1)-H.evalAt(cn,'neumann',1,t,1-dx))/dx<1e-3));
// 5 方程残差 u_t − k u_xx
let res=0;const k=0.7;for(const[x,t]of[[0.2,0.01],[0.47,0.03],[0.8,0.002]]){const h=1e-4,u=(x,t)=>H.evalAt(cm,'dirichlet',k,t,x);const ut=(u(x,t+h)-u(x,t-h))/(2*h),uxx=(u(x+h,t)-2*u(x,t)+u(x-h,t))/(h*h);res=Math.max(res,Math.abs(ut-k*uxx)/Math.max(1,Math.abs(ut)))}ok('级数解满足 u_t = k u_xx（差分残差）',res<1e-3,`相对残差 ${res.toExponential(2)}`);
// 6 Neumann 总热量守恒
let Q=0;const M=4000;for(let i=0;i<M;i++)Q+=H.evalAt(cn,'neumann',1,0.05,(i+.5)/M)/M;ok('绝热边界总热量守恒（三角初值 = 0.5）',Math.abs(Q-0.5)<1e-6,`∫u = ${Q.toFixed(8)}`);
// 7 热核闭式 vs 数值卷积
function conv(P,x,t,kk){let s=0;const h=0.0005,n=Math.round(24/h);for(let i=0;i<n;i++){const y=-12+(i+.5)*h;s+=H.kernel(x-y,t,kk)*H.linePhi(P,y)*h}return s}
let e7=0;for(const key of['box','two','gauss'])for(const[x,t]of[[0.3,0.05],[-1.2,0.4],[2,1]])e7=Math.max(e7,Math.abs(conv(H.linePresets[key],x,t,1)-H.lineU(H.linePresets[key],x,t,1)));ok('直线解闭式与 S*φ 数值卷积一致',e7<2e-6,`max 差 ${e7.toExponential(2)}`);
// 8 无穷传播速度
ok('方块初值在 t=0.01 时 x=5 处温度严格为正',H.lineU(H.linePresets.box,5,0.01,1)>0,`u = ${H.lineU(H.linePresets.box,5,0.01,1).toExponential(2)}`);
// 9 最大值原理：max_x u 随 t 不增
let prev=Infinity,mono=true;for(let lt=-3;lt<=0.7;lt+=0.05){const t=10**lt;let m=-1;for(let x=-5;x<=5;x+=0.01)m=Math.max(m,H.lineU(H.linePresets.two,x,t,1));if(m>prev+1e-9)mono=false;prev=m}ok('两块初值：max_x u 随时间单调不增',mono);
// 10 质量守恒
let mass=0;for(let x=-25;x<=25;x+=0.001)mass+=H.lineU(H.linePresets.box,x,0.5,1)*0.001;ok('直线上总热量守恒（方块 = 2）',Math.abs(mass-2)<1e-4,`∫u = ${mass.toFixed(6)}`);
// 11 von Neumann：单个离散模态一步后被乘以 G
const J=40;let e11=0;for(const r of[0.3,0.5,0.6])for(const m of[1,7,39]){const u=new Float64Array(J+1);for(let j=0;j<=J;j++)u[j]=Math.sin(m*PI*j/J);const v=H.fdStep(u,r),G=H.G(m*PI/J,r);for(let j=1;j<J;j++)e11=Math.max(e11,Math.abs(v[j]-G*u[j]))}ok('FTCS 对离散模态 sin(mπx_j) 的放大因子恰为 G(mπ/J)',e11<1e-12);
// 12 r ≤ ½ 稳定且收敛；r > ½ 爆炸
function run(r,steps){let u=new Float64Array(J+1);for(let j=0;j<=J;j++)u[j]=H.presets.tri(j/J);for(let i=0;i<steps;i++)u=H.fdStep(u,r);return u}
const r=0.4,dt=r/(J*J),n=200,us=run(r,n);let e12=0;for(let j=0;j<=J;j++)e12=Math.max(e12,Math.abs(us[j]-H.evalAt(ct,'dirichlet',1,n*dt,j/J)));ok('r = 0.4：差分解与级数解误差小',e12<5e-3,`max 误差 ${e12.toExponential(2)}`);
ok('r = 0.5：离散最大值原理（2000 步内 0 ≤ u ≤ 1）',run(0.5,2000).every(v=>v>=-1e-12&&v<=1+1e-12));
const blow=run(0.6,300);ok('r = 0.6：300 步内数值解爆炸（> 10³）',Math.max(...blow.map(Math.abs))>1e3,`max|u| = ${Math.max(...blow.map(Math.abs)).toExponential(2)}`);
console.log(`\n${pass} 项通过，${fail} 项失败`);process.exit(fail?1:0);
