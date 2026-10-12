// 数值核对：三个幂级数的系数之比趋于同一极限 l，且在 |z|<R 内通项趋于 0、|z|>R 外通项不趋于 0。
// 运行：node source/verify.cjs
const lnFact = n => { let s = 0; for (let k = 2; k <= n; k++) s += Math.log(k); return s; };
const presets = [
  { name: "a_n = 2^n",        l: 2,     R: 0.5,      r: n => 2,                     la: n => n * Math.LN2 },
  { name: "a_n = (n+1)/3^n",  l: 1 / 3, R: 3,        r: n => (n + 2) / (3 * (n + 1)), la: n => Math.log(n + 1) - n * Math.log(3) },
  { name: "a_n = 1/n!",       l: 0,     R: Infinity, r: n => 1 / (n + 1),           la: n => -lnFact(n) },
];
let pass = 0, fail = 0;
const check = (label, ok, detail) => { ok ? pass++ : fail++; console.log((ok ? "PASS " : "FAIL ") + label + (detail ? "  " + detail : "")); };

for (const p of presets) {
  // 系数之比：原级数 r(n)，逐项积分 r(n)(n+1)/(n+2)，逐项求导 r(n)(n+1)/n
  const ratios = [n => p.r(n), n => p.r(n) * (n + 1) / (n + 2), n => p.r(n) * (n + 1) / n];
  // r(n) 应与 la(n+1)-la(n) 一致
  let dev = 0;
  for (let n = 1; n <= 200; n++) dev = Math.max(dev, Math.abs(Math.log(p.r(n)) - (p.la(n + 1) - p.la(n))));
  check(`${p.name}: 比值公式与 ln|a_n| 一致`, dev < 1e-9, `最大偏差 ${dev.toExponential(2)}`);
  ratios.forEach((f, k) => {
    const N = 1e6, err = Math.abs(f(N) - p.l);
    check(`${p.name}: 第 ${k + 1} 个级数的系数之比在 n=1e6 处接近 l`, err < 1e-5, `|比值−l| = ${err.toExponential(2)}`);
  });
  // 通项的对数：ln|a_n z^n|、ln|a_n z^{n+1}/(n+1)|、ln|n a_n z^{n-1}|
  const terms = rho => [
    n => p.la(n) + n * Math.log(rho),
    n => p.la(n) - Math.log(n + 1) + (n + 1) * Math.log(rho),
    n => Math.log(n) + p.la(n) + (n - 1) * Math.log(rho),
  ];
  const inside = p.R === Infinity ? 50 : 0.9 * p.R;
  terms(inside).forEach((t, k) => check(`${p.name}: |z|=${inside} 在圆内，第 ${k + 1} 个级数的通项趋于 0`, t(4000) < -30 && t(4000) < t(2000)));
  if (p.R !== Infinity) {
    const outside = 1.1 * p.R;
    terms(outside).forEach((t, k) => check(`${p.name}: |z|=${outside.toFixed(2)} 在圆外，第 ${k + 1} 个级数的通项无界`, t(4000) > 30 && t(4000) > t(2000)));
  }
}
console.log(`\n${pass} 项通过，${fail} 项失败`);
process.exit(fail ? 1 : 0);
