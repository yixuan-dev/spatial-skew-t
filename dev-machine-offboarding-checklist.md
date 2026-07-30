# 開發機離職／交接 清除檢查清單

適用於 **Windows + Git + VS Code** 的開發機。離職、換機、把電腦交還 IT 或轉交同事之前照這份跑一遍。

> **核心觀念:Git 有三個彼此獨立的身分,清除時三個都要處理。**
>
> | 身分 | 由什麼決定 | 清除方式 |
> |---|---|---|
> | commit 作者 | `git config user.name/email` | 改 / 刪 gitconfig |
> | push 認證 | credential helper 快取的帳號 | 清認證管理員 |
> | GitHub 歸屬 | commit email 是否在該帳號驗證過 | (不需清,屬於遠端) |
>
> 只清其中一個 = 沒清乾淨。憑證是**還能繼續存取你帳號**的東西,不只是個資。

---

## 判斷:要做到什麼程度?

| 情境 | 需要做的範圍 |
|---|---|
| IT 整台重灌 / 系統重置 | 只需 **步驟 0**(確保工作不遺失),其餘可省 |
| 直接交接給下一個人 | **全部步驟**,尤其步驟 2 的憑證類 |
| 只是換自己的新機器 | 步驟 0 + 步驟 2,repo 可先備份再刪 |

---

## 步驟 0:先確保工作不會遺失(最重要,不可跳過)

刪任何東西之前先做完這一節。對**每個** repo 執行:

```powershell
$repo = "D:\Github\your-repo"       # 逐一替換

git -C $repo fetch --all
git -C $repo status -sb                              # 有沒有未 commit 的修改
git -C $repo log --branches --not --remotes --oneline # 所有分支上未推送的 commit
git -C $repo stash list                              # 忘在 stash 裡的東西
git -C $repo status --porcelain                      # 未追蹤檔案
```

四個指令**都要是空的 / 顯示同步**才算安全。

一次掃過所有 repo:

```powershell
Get-ChildItem D:\, "$env:USERPROFILE" -Filter .git -Recurse -Directory -Force -ErrorAction SilentlyContinue |
  ForEach-Object {
    $r = $_.Parent.FullName
    $unpushed = git -C $r log --branches --not --remotes --oneline 2>$null
    $dirty    = git -C $r status --porcelain 2>$null
    $stash    = git -C $r stash list 2>$null
    if ($unpushed -or $dirty -or $stash) {
      Write-Host "`n!! $r" -ForegroundColor Red
      if ($unpushed) { Write-Host "   未推送 commit: $($unpushed.Count)" }
      if ($dirty)    { Write-Host "   未 commit 變更: $($dirty.Count)" }
      if ($stash)    { Write-Host "   stash: $($stash.Count)" }
    } else { Write-Host "ok  $r" -ForegroundColor Green }
  }
```

> 遞迴掃描整顆磁碟可能要幾分鐘。

**別忘了不在 git 裡的東西:** 大型輸出檔、資料快取、`.gitignore` 掉的結果檔、桌面 / 下載資料夾的草稿、本機資料庫。這些 git 不會提醒你。

- [ ] 所有 repo 已推送、無 stash、無未 commit 變更
- [ ] git 以外的成果已另外備份

---

## 步驟 1:盤點(先看有什麼,再決定刪什麼)

每台機器留下的痕跡不同,先跑這段掃描,結果就是你的待刪清單。

```powershell
Write-Host "=== 1. 快取的憑證 ===" -ForegroundColor Cyan
cmdkey /list | Select-String -Pattern 'Target:|User:'

Write-Host "`n=== 2. Git 身分(三個層級) ===" -ForegroundColor Cyan
git config --show-origin --get-all user.name
git config --show-origin --get-all user.email

Write-Host "`n=== 3. 全域設定全文 ===" -ForegroundColor Cyan
git config --global --list

Write-Host "`n=== 4. 明文憑證檔(有的話是高風險) ===" -ForegroundColor Cyan
if (Test-Path "$env:USERPROFILE\.git-credentials") { Write-Host "!! 存在 ~/.git-credentials(明文 token)" -ForegroundColor Red } else { "ok - 無" }

Write-Host "`n=== 5. SSH 金鑰 ===" -ForegroundColor Cyan
if (Test-Path "$env:USERPROFILE\.ssh") { Get-ChildItem "$env:USERPROFILE\.ssh" | Select-Object Name } else { "ok - 無" }

Write-Host "`n=== 6. GPG 簽章金鑰 ===" -ForegroundColor Cyan
try { gpg --list-secret-keys 2>$null } catch { "ok - 無 gpg" }

Write-Host "`n=== 7. gh CLI 設定 ===" -ForegroundColor Cyan
if (Test-Path "$env:USERPROFILE\.config\gh") { Write-Host "!! 存在 ~/.config/gh(含 OAuth token)" -ForegroundColor Red } else { "ok - 無" }

Write-Host "`n=== 8. 環境變數裡的 token(只列名稱) ===" -ForegroundColor Cyan
Get-ChildItem Env: | Where-Object { $_.Name -match 'TOKEN|PAT|KEY|SECRET|PASS|CRED' } | Select-Object -ExpandProperty Name

Write-Host "`n=== 9. remote URL 裡有沒有內嵌 token ===" -ForegroundColor Cyan
Get-ChildItem D:\, "$env:USERPROFILE" -Filter .git -Recurse -Directory -Force -ErrorAction SilentlyContinue |
  ForEach-Object {
    $r = $_.Parent.FullName; $u = git -C $r remote get-url origin 2>$null
    if ($u -match '://[^/]*@') { Write-Host "!! $r -> URL 內嵌憑證" -ForegroundColor Red } else { "ok  $r -> $u" }
  }

Write-Host "`n=== 10. 工具的本機資料 ===" -ForegroundColor Cyan
@("$env:USERPROFILE\.claude", "$env:USERPROFILE\.claude.json", "$env:LOCALAPPDATA\Temp\claude",
  "$env:APPDATA\Code\User", "$env:USERPROFILE\.aws", "$env:USERPROFILE\.azure",
  "$env:USERPROFILE\.docker", "$env:USERPROFILE\.kube", "$env:USERPROFILE\.npmrc",
  "$env:USERPROFILE\.Renviron", "$env:USERPROFILE\.Rprofile") |
  ForEach-Object { if (Test-Path $_) { Write-Host "存在: $_" -ForegroundColor Yellow } }
```

把有輸出的項目抄下來,對照步驟 2 逐項處理。

---

## 步驟 2:清除

### 2.1 Git 認證

```powershell
cmdkey /delete:git:https://github.com
```

其他常見 target(依步驟 1 的結果調整):

```powershell
cmdkey /delete:git:https://gitlab.com
cmdkey /delete:git:https://dev.azure.com
cmdkey /delete:git:https://bitbucket.org
```

- [ ] Git 憑證已刪

### 2.2 Git 設定

```powershell
Remove-Item "$env:USERPROFILE\.gitconfig"
```

只想清身分、保留其他偏好設定:

```powershell
git config --global --unset user.name
git config --global --unset user.email
git config --global --unset-all credential.helper   # 視情況
```

> 系統層的 `C:\Program Files\Git\etc\gitconfig` 是 Git 安裝時的預設值,**不含個人資料,不用動**。

- [ ] 全域 gitconfig 已清

### 2.3 明文憑證與金鑰

```powershell
Remove-Item "$env:USERPROFILE\.git-credentials"  -ErrorAction SilentlyContinue
Remove-Item "$env:USERPROFILE\.ssh"    -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item "$env:USERPROFILE\.config\gh" -Recurse -Force -ErrorAction SilentlyContinue
```

有 GPG 簽章金鑰的話:

```powershell
gpg --list-secret-keys --keyid-format=long      # 找出 key id
gpg --delete-secret-and-public-key <KEY_ID>
```

- [ ] SSH / GPG / gh token 已清

### 2.4 Repo 本身

**確認步驟 0 全部通過後**才執行(不可逆):

```powershell
Remove-Item -Recurse -Force D:\Github\your-repo
```

repo 的 local config(`.git\config`,含 per-repo 身分與 remote)會跟著一起消失,不用另外處理。

- [ ] repo 已刪除或已備份帶走

### 2.5 VS Code

**用 UI 登出,不要只刪檔案** —— VS Code 的登入 token 存在自己的 secret storage,`cmdkey` 清不到:

- 左下角**帳號圖示 → GitHub → 登出**
- 其他有登入的擴充套件(Copilot、Claude Code、Docker、雲端服務)一併登出

要整個重置設定與擴充套件資料:

```powershell
Remove-Item -Recurse -Force "$env:APPDATA\Code\User"
Remove-Item -Recurse -Force "$env:USERPROFILE\.vscode\extensions"
```

- [ ] VS Code 及擴充套件已登出

### 2.6 AI 工具的本機資料

```powershell
Remove-Item -Recurse -Force "$env:USERPROFILE\.claude"
Remove-Item -Force        "$env:USERPROFILE\.claude.json" -ErrorAction SilentlyContinue
Remove-Item -Recurse -Force "$env:LOCALAPPDATA\Temp\claude"
```

> `.claude\projects`、`file-history`、`shell-snapshots` 存的是**對話記錄與程式碼快照**,內容通常比預期多,別漏。

- [ ] AI 工具資料已清

### 2.7 其他開發者憑證(有裝才需要)

```powershell
Remove-Item -Recurse -Force "$env:USERPROFILE\.aws"     -ErrorAction SilentlyContinue
Remove-Item -Recurse -Force "$env:USERPROFILE\.azure"   -ErrorAction SilentlyContinue
Remove-Item -Recurse -Force "$env:USERPROFILE\.kube"    -ErrorAction SilentlyContinue
Remove-Item -Recurse -Force "$env:USERPROFILE\.docker"  -ErrorAction SilentlyContinue
Remove-Item -Force "$env:USERPROFILE\.npmrc"            -ErrorAction SilentlyContinue
Remove-Item -Force "$env:USERPROFILE\.pypirc"           -ErrorAction SilentlyContinue
Remove-Item -Force "$env:USERPROFILE\.Renviron"         -ErrorAction SilentlyContinue
```

`.npmrc` / `.pypirc` / `.Renviron` 常藏 registry token 或 `GITHUB_PAT`,容易被忽略。

- [ ] 其他 CLI 憑證已清

### 2.8 雲端硬碟 / 作業系統帳號

**用各 App 自己的登出功能,不要直接刪認證管理員項目** —— 直接刪會讓 App 卡在半登入狀態,而且本機同步的檔案快取還留在硬碟上。

| 服務 | 做法 |
|---|---|
| Google Drive | 偏好設定 → 中斷帳戶連結 |
| OneDrive | 設定 → 帳戶 → 解除此電腦的連結 |
| Microsoft 帳號 | 設定 → 帳戶 → 中斷連結 |

- [ ] 雲端硬碟已解除連結且本機快取已清

### 2.9 瀏覽器

- 登出所有工作相關網站(GitHub 尤其重要)
- 清除儲存的密碼與 cookie
- 移除同步中的設定檔

> **踩過的坑:** 瀏覽器留著登入 session 的話,之後 Git 走 OAuth 流程會**默默沿用**那個帳號,等於留了一把鑰匙。

- [ ] 瀏覽器已登出並清除

---

## 步驟 3:驗證

```powershell
Write-Host "=== 憑證(應無輸出) ===" -ForegroundColor Cyan
cmdkey /list | Select-String github, gitlab, azure

Write-Host "`n=== Git 身分(應無 user.*) ===" -ForegroundColor Cyan
git config --global --list

Write-Host "`n=== 檔案(應全為 False) ===" -ForegroundColor Cyan
@("$env:USERPROFILE\.gitconfig", "$env:USERPROFILE\.ssh", "$env:USERPROFILE\.git-credentials",
  "$env:USERPROFILE\.claude", "$env:USERPROFILE\.config\gh") |
  ForEach-Object { "{0,-50} {1}" -f $_, (Test-Path $_) }
```

最後在遠端側收尾(從別台機器或手機做):

- [ ] GitHub → Settings → **Sessions** 撤銷這台機器的 session
- [ ] GitHub → Settings → **Applications / Authorized OAuth Apps** 撤銷 Git Credential Manager 授權
- [ ] GitHub → Settings → **Personal access tokens** 撤銷這台機器用過的 token
- [ ] GitHub → Settings → **SSH and GPG keys** 刪掉這台機器的公鑰

> 這一步最容易漏。本機刪掉憑證只是**你這台**登不進去;遠端的 token 與授權若沒撤銷,那把鑰匙依然有效。

---

## 附錄 A:常見陷阱

| 陷阱 | 說明 |
|---|---|
| 只改了 `git config` 就以為清乾淨 | 作者身分與 push 認證是**兩套獨立系統**,改 config 動不到憑證 |
| `Permission to X denied to Y` 誤判 | **Y 是認證身分**,不是 commit 作者。看到這句直接查 credential helper |
| VS Code 帳號沒登出 | token 在 VS Code 自己的 secret storage,`cmdkey` 掃不到也刪不掉 |
| 瀏覽器 session 沒清 | OAuth 會默默沿用,等於憑證沒刪 |
| 遠端授權沒撤銷 | 本機清乾淨 ≠ 帳號安全,token 在伺服器端依然有效 |
| 忘記 `.gitignore` 掉的成果 | 大型輸出、資料快取 git 不會提醒,一刪就沒了 |
| 用 `--global` 設身分導致混用 | 多帳號時應該用 per-repo 的 local config |

## 附錄 B:macOS / Linux 對照

| 項目 | Windows | macOS | Linux |
|---|---|---|---|
| 憑證儲存 | 認證管理員 (`cmdkey`) | 鑰匙圈 (Keychain Access) | libsecret / `~/.git-credentials` |
| 刪 Git 憑證 | `cmdkey /delete:git:https://github.com` | 鑰匙圈搜尋 `github.com` → 刪除 | `secret-tool clear service github.com` |
| 全域設定 | `%USERPROFILE%\.gitconfig` | `~/.gitconfig` | `~/.gitconfig` |
| SSH 金鑰 | `%USERPROFILE%\.ssh` | `~/.ssh` | `~/.ssh` |
| 通用登出 | — | `git credential-osxkeychain erase` | `git credential-cache exit` |

跨平台通用(任何 OS 都可用):

```bash
printf "protocol=https\nhost=github.com\n\n" | git credential reject
```

## 附錄 C:新機器上線時(反向清單)

換到新機器時,**第一個 commit 之前**先做這三件事,就不會重蹈覆轍:

```bash
git config --global user.name  "你的名字"
git config --global user.email "你的@email"
git config --global --get user.email     # 確認
```

- 多帳號並存時,用 per-repo 的 local config(不加 `--global`)
- email 必須是**在該 GitHub 帳號下驗證過**的,否則 commit 不會連結到你的帳號、也不計入 contribution
- 考慮改用 SSH:一組 key 綁一個帳號,不受瀏覽器 session 干擾

---

*用途:個人開發機交接清除。使用前先跑步驟 0,刪除操作不可逆。*
