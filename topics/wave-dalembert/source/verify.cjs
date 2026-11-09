// 运行：node source/verify.cjs —— 从 index.html 抽取 <script id="core"> 的纯数学内核并做独立数值检查
const fs=require('fs'),path=require('path');
const html=fs.readFileSync(path.join(__dirname,'..','index.html'),'utf8');
const W=new Function(html.match(/<script id="core">([\s\S]*?)<\/script>/)[1]+'\nreturn WaveCore;')();
let pass=0,fail=0;const ok=(name,cond,info='')=>{if(cond){pass++;console.log('✓',name,info)}else{fail++;console.log('✗',name,info)}};
const G=W.phiPresets.gauss,P=W.psiPresets.gauss,c=1.3,h=1e-3;
// 1 初值条件
let e1=0,e2=0;for(let x=-3;x<=3;x+=0.37){e1=Math.max(e1,Math.abs(W.u(G,P,c,x,0)-G.f(x)));const ut=(W.u(G,P,c,x,h)-W.u(G,P,c,x,-h))/(2*h);e2=Math.max(e2,Math.abs(ut-P.f(x)))}
ok('u(x,0) = φ(x)',e1<1e-12);ok('u_t(x,0) = ψ(x)（中心差分）',e2<1e-5,`误差 ${e2.toExponential(2)}`);
// 2 方程残差
let res=0;const hh=2e-3;for(const[x,t]of[[0.3,0.5],[-1.1,1.2],[2,0.9],[0,2.5]]){const u=(x,t)=>W.u(G,P,c,x,t);const utt=(u(x,t+hh)-2*u(x,t)+u(x,t-hh))/hh**2,uxx=(u(x+hh,t)-2*u(x,t)+u(x-hh,t))/hh**2;res=Math.max(res,Math.abs(utt-c*c*uxx))}
ok('d\'Alembert 解满足 u_tt = c²u_xx（差分残差）',res<1e-4,`残差 ${res.toExponential(2)}`);
// 3 原函数正确：F′ = f
let e3=0;for(const k of['gauss','box'])for(let x=-2;x<=2;x+=0.173){const p=W.psiPresets[k];e3=Math.max(e3,Math.abs((p.F(x+1e-5)-p.F(x-1e-5))/2e-5-p.f(x)))}ok('ψ 的原函数 F 满足 F′ = ψ',e3<1e-5);
// 4 方块冲量留下的平台高度 = ∫ψ/(2c)
ok('方块冲量：t 足够大时 u(0,t) = 1.4/(2c)',Math.abs(W.u(W.phiPresets.none,W.psiPresets.box,c,0,3)-1.4/(2*c))<1e-12);
// 5 有限弦：边界条件、周期、能量、反射符号
const Lw=4;let bcF=0,bcN=0,per=0,dE=0;
for(const kind of['disp','vel'])for(const bc of['fixed','free']){const F=W.finite({Lw,p:1.2,w:0.3,A:1,kind,bc,c:1});const E0=F.energy(0,3000).tot;
  for(let t=0;t<=16;t+=0.731){if(bc==='fixed')bcF=Math.max(bcF,Math.abs(F.U(0,t)),Math.abs(F.U(Lw,t)));else bcN=Math.max(bcN,Math.abs(F.Ux(0,t)),Math.abs(F.Ux(Lw,t)));dE=Math.max(dE,Math.abs(F.energy(t,3000).tot-E0)/E0)}
  if(bc==='fixed')for(let x=0.1;x<Lw;x+=0.37)per=Math.max(per,Math.abs(F.U(x,1.7+8)-F.U(x,1.7)))}
ok('固定端：u(0,t) = u(L,t) = 0',bcF<1e-12,`max ${bcF.toExponential(2)}`);ok('自由端：uₓ(0,t) = uₓ(L,t) = 0',bcN<1e-10,`max ${bcN.toExponential(2)}`);
ok('两端固定：解以 2L/c 为周期',per<1e-12);ok('四种组合下能量守恒（相对漂移 < 1e-9）',dE<1e-9,`max ${dE.toExponential(2)}`);
const Ff=W.finite({Lw,p:1.2,w:0.3,A:1,kind:'disp',bc:'fixed',c:1}),Fn=W.finite({Lw,p:1.2,w:0.3,A:1,kind:'disp',bc:'free',c:1});
ok('固定端反射后翻转：u(p, 2p) ≈ −½',Math.abs(Ff.U(1.2,2.4)+0.5)<1e-3,`u = ${Ff.U(1.2,2.4).toFixed(5)}`);ok('自由端反射后不翻转：u(p, 2p) ≈ +½',Math.abs(Fn.U(1.2,2.4)-0.5)<1e-3,`u = ${Fn.U(1.2,2.4).toFixed(5)}`);
// 6 解析导数与数值导数一致
let e6=0;for(const F of[Ff,Fn])for(const[x,t]of[[0.7,1.1],[2.3,3.3],[3.9,5]]){e6=Math.max(e6,Math.abs(F.Ut(x,t)-(F.U(x,t+1e-5)-F.U(x,t-1e-5))/2e-5),Math.abs(F.Ux(x,t)-(F.U(x+1e-5,t)-F.U(x-1e-5,t))/2e-5))}ok('u_t、u_x 解析式与数值导数一致',e6<1e-5,`max ${e6.toExponential(2)}`);
console.log(`\n${pass} 项通过，${fail} 项失败`);process.exit(fail?1:0);
