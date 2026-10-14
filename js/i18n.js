/* i18n.js — trilingual UI strings (EN / zh-Hans / zh-Hant).
 * Usage: window.setLang(code); window.t(key) returns the current-locale string.
 */
(function (root) {
  "use strict";

  var DICT = {
    en: {
      brand: "BioAI Lab",
      product: "LoCRANE · Aging Clock",
      hero: "Drop in a transcriptomic dataset and get a longevity-certified aging clock — computed entirely in your browser.",
      heroSub: "No upload, no server. Your expression matrix + phenotype never leave this page.",
      dropTitle: "Drag your CSVs here",
      dropOr: "or",
      dropSub: "Expression matrix (genes × samples) + phenotype (sample_id, age, sex). Comma- or tab-separated.",
      btnDemo: "Load built-in demo",
      btnDownload: "Download demo CSVs",
      dropHintExpr: "expression matrix (genes × samples)",
      dropHintPhen: "phenotype (sample_id, age, sex)",
      computing: "Computing in your browser…",
      resultsTitle: "Results",
      resultsSub: "Linear core · 5-fold cross-validated · computed locally",
      kpiMae: "Ensemble MAE",
      kpiR2: "R²",
      kpiSpearman: "Spearman ρ",
      kpiSamples: "Samples",
      kpiGenes: "Genes",
      kpiSex: "Sex prior value",
      kpiGate: "90+ reads younger",
      viewsTitle: "Model views · 5-fold CV MAE (lower = better)",
      viewAll: "Ridge · all genes",
      view40: "EN · top 40",
      viewBig: "EN · top 150",
      viewEns: "Ensemble",
      scatterTitle: "Predicted biological age vs chronological age",
      scatterX: "Chronological age",
      scatterY: "Predicted age",
      bandTitle: "Per-band MAE · Y / M / L",
      bandY: "Y · ≤50",
      bandM: "M · 51–89",
      bandL: "L · ≥90",
      leversTitle: "Top causal levers · |β|",
      leversSub: "Genes with the largest standardized weights in the all-gene ridge clock",
      gateTitle: "90+ longevity gate",
      gateAuc: "Gate AUC",
      gateMirror: "90+ tail is picked out by its young read-out (mirror AUC {m})",
      gateNote: "The 90+ tail reads {v} years younger than its chronological age: reaching 90 (the gate) decouples from how old the blood looks (the rate).",
      sampleTable: "Per-sample read-out",
      colSample: "Sample",
      colChron: "Chron age",
      colPred: "Pred age",
      colResid: "Residual",
      colBand: "Band",
      disclaimer: "Demo-scale linear core (ridge + elastic-net feature budgets, 5-fold CV), computed locally in your browser. The full deep LoCRANE ensemble runs via code/run_demo.py in the repository.",
      footerPaper: "LoCRANE · Engineering Applications of Artificial Intelligence",
      footerRepo: "Repository",
      langLabel: "Language",
      errParse: "Couldn't parse that. Use a genes×samples expression CSV and a sample_id,age,sex phenotype CSV.",
      errMatch: "Fewer than 8 expression samples matched the phenotype IDs.",
      errUnknown: "Something went wrong computing the clock."
    },
    "zh-Hans": {
      brand: "BioAI Lab",
      product: "LoCRANE · 老化时钟",
      hero: "拖入一个转录组数据集，即可得到一枚「长寿认证」的老化时钟——全部在你的浏览器里完成计算。",
      heroSub: "不上传、不走服务器。你的表达矩阵与表型不会离开当前页面。",
      dropTitle: "把 CSV 拖到这里",
      dropOr: "或",
      dropSub: "表达矩阵（基因 × 样本）＋ 表型（sample_id, age, sex）。逗号或制表符分隔。",
      btnDemo: "载入内置示例",
      btnDownload: "下载示例 CSV",
      dropHintExpr: "表达矩阵（基因 × 样本）",
      dropHintPhen: "表型（sample_id, age, sex）",
      computing: "正在你的浏览器里计算…",
      resultsTitle: "结果",
      resultsSub: "线性核 · 五折交叉验证 · 本地计算",
      kpiMae: "集成 MAE",
      kpiR2: "R²",
      kpiSpearman: "Spearman ρ",
      kpiSamples: "样本数",
      kpiGenes: "基因数",
      kpiSex: "性别先验价值",
      kpiGate: "90+ 读得更年轻",
      viewsTitle: "各模型视图 · 五折 CV MAE（越低越好）",
      viewAll: "岭回归 · 全基因",
      view40: "弹网 · 前 40",
      viewBig: "弹网 · 前 150",
      viewEns: "集成",
      scatterTitle: "预测生物学年龄 vs 日历年龄",
      scatterX: "日历年龄",
      scatterY: "预测年龄",
      bandTitle: "分带 MAE · Y / M / L",
      bandY: "Y · ≤50",
      bandM: "M · 51–89",
      bandL: "L · ≥90",
      leversTitle: "Top 因果杠杆基因 · |β|",
      leversSub: "在全基因岭回归钟里标准化权重最大的基因",
      gateTitle: "90+ 长寿之门",
      gateAuc: "门 AUC",
      gateMirror: "90+ 尾端由『年轻读数』分出（镜像 AUC {m}）",
      gateNote: "90+ 这一尾端读出的年龄比其日历年龄年轻 {v} 岁：够到 90（门）与血液看起来多老（速率）彼此解耦。",
      sampleTable: "逐样本回答",
      colSample: "样本",
      colChron: "日历年龄",
      colPred: "预测年龄",
      colResid: "残差",
      colBand: "带",
      disclaimer: "演示规模的线性核（岭回归＋弹性网特征预算，五折交叉验证），完全在你的浏览器本地计算。完整的 LoCRANE 深度集成模型请在仓库 code/run_demo.py 中运行。",
      footerPaper: "LoCRANE · Engineering Applications of Artificial Intelligence",
      footerRepo: "代码仓库",
      langLabel: "语言",
      errParse: "无法解析。请使用 基因×样本 的表达 CSV 与 sample_id,age,sex 的表型 CSV。",
      errMatch: "表达样本中匹配到表型 ID 的不足 8 个。",
      errUnknown: "计算时钟时出错了。"
    },
    "zh-Hant": {
      brand: "BioAI Lab",
      product: "LoCRANE · 老化時鐘",
      hero: "拖入一個轉錄組資料集，即可得到一枚「長壽認證」的老化時鐘——全部在你的瀏覽器裡完成計算。",
      heroSub: "不上傳、不走伺服器。你的表達矩陣與表現型不會離開當前頁面。",
      dropTitle: "把 CSV 拖到這裡",
      dropOr: "或",
      dropSub: "表達矩陣（基因 × 樣本）＋ 表現型（sample_id, age, sex）。逗號或直欄分隔。",
      btnDemo: "載入內建範例",
      btnDownload: "下載範例 CSV",
      dropHintExpr: "表達矩陣（基因 × 樣本）",
      dropHintPhen: "表現型（sample_id, age, sex）",
      computing: "正在你的瀏覽器裡計算…",
      resultsTitle: "結果",
      resultsSub: "線性核 · 五折交叉驗證 · 本機計算",
      kpiMae: "集成 MAE",
      kpiR2: "R²",
      kpiSpearman: "Spearman ρ",
      kpiSamples: "樣本數",
      kpiGenes: "基因數",
      kpiSex: "性別先驗價值",
      kpiGate: "90+ 讀得更年輕",
      viewsTitle: "各模型視圖 · 五折 CV MAE（越低越好）",
      viewAll: "稜迴歸 · 全基因",
      view40: "彈網 · 前 40",
      viewBig: "彈網 · 前 150",
      viewEns: "集成",
      scatterTitle: "預測生物學年齡 vs 日曆年齡",
      scatterX: "日曆年齡",
      scatterY: "預測年齡",
      bandTitle: "分帶 MAE · Y / M / L",
      bandY: "Y · ≤50",
      bandM: "M · 51–89",
      bandL: "L · ≥90",
      leversTitle: "Top 因果槓桿基因 · |β|",
      leversSub: "在全基因稜迴歸鐘中標準化權重最大的基因",
      gateTitle: "90+ 長壽之門",
      gateAuc: "門 AUC",
      gateMirror: "90+ 尾端由『年輕讀數』分出（鏡像 AUC {m}）",
      gateNote: "90+ 這個尾端讀出的年齡比其日曆年齡年輕 {v} 歲：夠到 90（門）與血液看起來多老（速率）彼此解耦。",
      sampleTable: "逐樣本答讀",
      colSample: "樣本",
      colChron: "日曆年齡",
      colPred: "預測年齡",
      colResid: "殘差",
      colBand: "帶",
      disclaimer: "示範規模的線性核（稜迴歸＋彈性網特徵預算，五折交叉驗證），完全在你的瀏覽器本機計算。完整的 LoCRANE 深度集成模型請在倉庫 code/run_demo.py 中執行。",
      footerPaper: "LoCRANE · Engineering Applications of Artificial Intelligence",
      footerRepo: "程式碼倉庫",
      langLabel: "語言",
      errParse: "無法解析。請使用 基因×樣本 的表達 CSV 與 sample_id,age,sex 的表現型 CSV。",
      errMatch: "表達樣本中對到表現型 ID 的不足 8 個。",
      errUnknown: "計算時鐘時出錯。"
    }
  };

  var current = "en";
  function detect() {
    var s = root.localStorage ? root.localStorage.getItem("bioai-lang") : null;
    if (s && DICT[s]) return s;
    return "en";
  }
  current = detect();

  root.I18N = DICT;
  root.getLang = function () { return current; };
  root.setLang = function (code) {
    if (!DICT[code]) code = "en";
    current = code;
    if (root.localStorage) { try { root.localStorage.setItem("bioai-lang", code); } catch (e) {} }
    applyDom();
  };
  root.t = function (key, vars) {
    var s = (DICT[current] && DICT[current][key]) !== undefined ? DICT[current][key]
      : (DICT.en[key] !== undefined ? DICT.en[key] : key);
    if (vars) Object.keys(vars).forEach(function (k) { s = s.replace(new RegExp("\\{" + k + "\\}", "g"), vars[k]); });
    return s;
  };

  function applyDom() {
    var byId = {};
    document.querySelectorAll("[data-i18n]").forEach(function (el) {
      byId[el.getAttribute("data-i18n")] = el;
      el.textContent = t(el.getAttribute("data-i18n"));
    });
    if (window.__onLang) window.__onLang(current, byId);
  }
  root.__applyDom = applyDom;
  root.__onLang = null;

  root.__langCodes = ["en", "zh-Hans", "zh-Hant"];
  root.__langLabel = { "en": "EN", "zh-Hans": "中文", "zh-Hant": "粵語" };
})(typeof window !== "undefined" ? window : this);
