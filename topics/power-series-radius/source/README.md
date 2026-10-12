# 源码与验证

专题页 [`../index.html`](../index.html) 本身就是可编辑的源码：样式、脚本和正文都在这一个文件里，没有构建步骤。

- `verify.cjs`：用页面里同样的系数公式做数值核对。

```bash
node source/verify.cjs
```

预览图用无头 Chrome 生成（在专题目录下运行）：

```bash
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless=new --hide-scrollbars --window-size=1200,675 --virtual-time-budget=8000 --screenshot="$PWD/preview.png" "file://$PWD/index.html"
```
