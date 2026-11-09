# 热方程专题 · 源码与验证

`../index.html` 本身就是可编辑源码（单文件，无构建步骤）。其中 `<script id="core">` 是不依赖 DOM 的数学内核 `HeatCore`，页面交互代码在其后的 `<script>` 中。

```bash
node source/verify.cjs   # 从 index.html 抽取 HeatCore 并做 15 项数值检查
```
