#if defined(__ELF__) && defined(__linux__)
.section .note.GNU-stack,"",%progbits
#endif

#if defined(__ELF__) && defined(__CET__) && defined(__has_include)
#if __has_include(<cet.h>)
#include <cet.h>
#endif
#endif

#if !defined(_CET_ENDBR)
#define _CET_ENDBR
#endif

.intel_syntax noprefix
.global _blake3_hash_many_avx2
.global blake3_hash_many_avx2
#ifdef __APPLE__
.text
#else
.section .text
#endif
        .p2align  6
_blake3_hash_many_avx2:
blake3_hash_many_avx2:
        _CET_ENDBR
        push    r15
        push    r14
        push    r13
        push    r12
        push    rbx
        push    rbp
        mov     rbp, rsp
        sub     rsp, 680
        and     rsp, 0xFFFFFFFFFFFFFFC0
        neg     r9d
        vmovd   xmm0, r9d
        vpbroadcastd ymm0, xmm0
        vmovdqa ymmword ptr [rsp+0x280], ymm0
        vpand   ymm1, ymm0, ymmword ptr [ADD0+rip]
        vpand   ymm2, ymm0, ymmword ptr [ADD1+rip]
        vmovdqa ymmword ptr [rsp+0x220], ymm2
        vmovd   xmm2, r8d
        vpbroadcastd ymm2, xmm2
        vpaddd  ymm2, ymm2, ymm1
        vmovdqa ymmword ptr [rsp+0x240], ymm2
        vpxor   ymm1, ymm1, ymmword ptr [CMP_MSB_MASK+rip]
        vpxor   ymm2, ymm2, ymmword ptr [CMP_MSB_MASK+rip]
        vpcmpgtd ymm2, ymm1, ymm2
        shr     r8, 32
        vmovd   xmm3, r8d
        vpbroadcastd ymm3, xmm3
        vpsubd  ymm3, ymm3, ymm2
        vmovdqa ymmword ptr [rsp+0x260], ymm3
        shl     rdx, 6
        mov     qword ptr [rsp+0x2A0], rdx
        cmp     rsi, 8
        jc      3f
2:
        vpbroadcastd ymm0, dword ptr [rcx]
        vpbroadcastd ymm1, dword ptr [rcx+0x4]
        vpbroadcastd ymm2, dword ptr [rcx+0x8]
        vpbroadcastd ymm3, dword ptr [rcx+0xC]
        vpbroadcastd ymm4, dword ptr [rcx+0x10]
        vpbroadcastd ymm5, dword ptr [rcx+0x14]
        vpbroadcastd ymm6, dword ptr [rcx+0x18]
        vpbroadcastd ymm7, dword ptr [rcx+0x1C]
        mov     r8, qword ptr [rdi]
        mov     r9, qword ptr [rdi+0x8]
        mov     r10, qword ptr [rdi+0x10]
        mov     r11, qword ptr [rdi+0x18]
        mov     r12, qword ptr [rdi+0x20]
        mov     r13, qword ptr [rdi+0x28]
        mov     r14, qword ptr [rdi+0x30]
        mov     r15, qword ptr [rdi+0x38]
        movzx   eax, byte ptr [rbp+0x38]
        movzx   ebx, byte ptr [rbp+0x40]
        or      eax, ebx
        xor     edx, edx
.p2align  5
9:
        movzx   ebx, byte ptr [rbp+0x48]
        or      ebx, eax
        add     rdx, 64
        cmp     rdx, qword ptr [rsp+0x2A0]
        cmove   eax, ebx
        mov     dword ptr [rsp+0x200], eax
        vmovups xmm8, xmmword ptr [r8+rdx-0x40]
        vinsertf128 ymm8, ymm8, xmmword ptr [r12+rdx-0x40], 0x01
        vmovups xmm9, xmmword ptr [r9+rdx-0x40]
        vinsertf128 ymm9, ymm9, xmmword ptr [r13+rdx-0x40], 0x01
        vunpcklpd ymm12, ymm8, ymm9
        vunpckhpd ymm13, ymm8, ymm9
        vmovups xmm10, xmmword ptr [r10+rdx-0x40]
        vinsertf128 ymm10, ymm10, xmmword ptr [r14+rdx-0x40], 0x01
        vmovups xmm11, xmmword ptr [r11+rdx-0x40]
        vinsertf128 ymm11, ymm11, xmmword ptr [r15+rdx-0x40], 0x01
        vunpcklpd ymm14, ymm10, ymm11
        vunpckhpd ymm15, ymm10, ymm11
        vshufps ymm8, ymm12, ymm14, 136
        vmovaps ymmword ptr [rsp], ymm8
        vshufps ymm9, ymm12, ymm14, 221
        vmovaps ymmword ptr [rsp+0x20], ymm9
        vshufps ymm10, ymm13, ymm15, 136
        vmovaps ymmword ptr [rsp+0x40], ymm10
        vshufps ymm11, ymm13, ymm15, 221
        vmovaps ymmword ptr [rsp+0x60], ymm11
        vmovups xmm8, xmmword ptr [r8+rdx-0x30]
        vinsertf128 ymm8, ymm8, xmmword ptr [r12+rdx-0x30], 0x01
        vmovups xmm9, xmmword ptr [r9+rdx-0x30]
        vinsertf128 ymm9, ymm9, xmmword ptr [r13+rdx-0x30], 0x01
        vunpcklpd ymm12, ymm8, ymm9
        vunpckhpd ymm13, ymm8, ymm9
        vmovups xmm10, xmmword ptr [r10+rdx-0x30]
        vinsertf128 ymm10, ymm10, xmmword ptr [r14+rdx-0x30], 0x01
        vmovups xmm11, xmmword ptr [r11+rdx-0x30]
        vinsertf128 ymm11, ymm11, xmmword ptr [r15+rdx-0x30], 0x01
        vunpcklpd ymm14, ymm10, ymm11
        vunpckhpd ymm15, ymm10, ymm11
        vshufps ymm8, ymm12, ymm14, 136
        vmovaps ymmword ptr [rsp+0x80], ymm8
        vshufps ymm9, ymm12, ymm14, 221
        vmovaps ymmword ptr [rsp+0xA0], ymm9
        vshufps ymm10, ymm13, ymm15, 136
        vmovaps ymmword ptr [rsp+0xC0], ymm10
        vshufps ymm11, ymm13, ymm15, 221
        vmovaps ymmword ptr [rsp+0xE0], ymm11
        vmovups xmm8, xmmword ptr [r8+rdx-0x20]
        vinsertf128 ymm8, ymm8, xmmword ptr [r12+rdx-0x20], 0x01
        vmovups xmm9, xmmword ptr [r9+rdx-0x20]
        vinsertf128 ymm9, ymm9, xmmword ptr [r13+rdx-0x20], 0x01
        vunpcklpd ymm12, ymm8, ymm9
        vunpckhpd ymm13, ymm8, ymm9
        vmovups xmm10, xmmword ptr [r10+rdx-0x20]
        vinsertf128 ymm10, ymm10, xmmword ptr [r14+rdx-0x20], 0x01
        vmovups xmm11, xmmword ptr [r11+rdx-0x20]
        vinsertf128 ymm11, ymm11, xmmword ptr [r15+rdx-0x20], 0x01
        vunpcklpd ymm14, ymm10, ymm11
        vunpckhpd ymm15, ymm10, ymm11
        vshufps ymm8, ymm12, ymm14, 136
        vmovaps ymmword ptr [rsp+0x100], ymm8
        vshufps ymm9, ymm12, ymm14, 221
        vmovaps ymmword ptr [rsp+0x120], ymm9
        vshufps ymm10, ymm13, ymm15, 136
        vmovaps ymmword ptr [rsp+0x140], ymm10
        vshufps ymm11, ymm13, ymm15, 221
        vmovaps ymmword ptr [rsp+0x160], ymm11
        vmovups xmm8, xmmword ptr [r8+rdx-0x10]
        vinsertf128 ymm8, ymm8, xmmword ptr [r12+rdx-0x10], 0x01
        vmovups xmm9, xmmword ptr [r9+rdx-0x10]
        vinsertf128 ymm9, ymm9, xmmword ptr [r13+rdx-0x10], 0x01
        vunpcklpd ymm12, ymm8, ymm9
        vunpckhpd ymm13, ymm8, ymm9
        vmovups xmm10, xmmword ptr [r10+rdx-0x10]
        vinsertf128 ymm10, ymm10, xmmword ptr [r14+rdx-0x10], 0x01
        vmovups xmm11, xmmword ptr [r11+rdx-0x10]
        vinsertf128 ymm11, ymm11, xmmword ptr [r15+rdx-0x10], 0x01
        vunpcklpd ymm14, ymm10, ymm11
        vunpckhpd ymm15, ymm10, ymm11
        vshufps ymm8, ymm12, ymm14, 136
        vmovaps ymmword ptr [rsp+0x180], ymm8
        vshufps ymm9, ymm12, ymm14, 221
        vmovaps ymmword ptr [rsp+0x1A0], ymm9
        vshufps ymm10, ymm13, ymm15, 136
        vmovaps ymmword ptr [rsp+0x1C0], ymm10
        vshufps ymm11, ymm13, ymm15, 221
        vmovaps ymmword ptr [rsp+0x1E0], ymm11
        vpbroadcastd ymm15, dword ptr [rsp+0x200]
        prefetcht0 [r8+rdx+0x80]
        prefetcht0 [r12+rdx+0x80]
        prefetcht0 [r9+rdx+0x80]
        prefetcht0 [r13+rdx+0x80]
        prefetcht0 [r10+rdx+0x80]
        prefetcht0 [r14+rdx+0x80]
        prefetcht0 [r11+rdx+0x80]
        prefetcht0 [r15+rdx+0x80]
        vpaddd  ymm0, ymm0, ymmword ptr [rsp]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x40]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x80]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0xC0]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm0, ymmword ptr [rsp+0x240]
        vpxor   ymm13, ymm1, ymmword ptr [rsp+0x260]
        vpxor   ymm14, ymm2, ymmword ptr [BLAKE3_BLOCK_LEN+rip]
        vpxor   ymm15, ymm3, ymm15
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [BLAKE3_IV_0+rip]
        vpaddd  ymm9, ymm13, ymmword ptr [BLAKE3_IV_1+rip]
        vpaddd  ymm10, ymm14, ymmword ptr [BLAKE3_IV_2+rip]
        vpaddd  ymm11, ymm15, ymmword ptr [BLAKE3_IV_3+rip]
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x20]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x60]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0xA0]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0xE0]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x100]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x140]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x180]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x1C0]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x120]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x160]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x1A0]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x1E0]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x40]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x60]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0xE0]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x80]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0xC0]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x140]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x1A0]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x20]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x180]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x120]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x1E0]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x160]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0xA0]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x1C0]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x100]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x60]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x140]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x1A0]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0xE0]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x80]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x180]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x40]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x1C0]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0xC0]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x120]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x160]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x100]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0xA0]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x1E0]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x20]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x140]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x180]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x1C0]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x1A0]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0xE0]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x120]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x60]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x1E0]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x80]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x160]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0xA0]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x20]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x40]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x100]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0xC0]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x180]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x120]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x1E0]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x1C0]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x1A0]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x160]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x140]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x100]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0xE0]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0xA0]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0xC0]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x40]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x60]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x20]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x80]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x120]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x160]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x100]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x1E0]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x1C0]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0xA0]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x180]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x20]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x1A0]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x40]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x80]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x60]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x140]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0xC0]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0xE0]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x160]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0xA0]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x20]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x100]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x1E0]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x120]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0xC0]
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm1, ymm1, ymm5
        vpaddd  ymm2, ymm2, ymm6
        vpaddd  ymm3, ymm3, ymm7
        vpxor   ymm12, ymm12, ymm0
        vpxor   ymm13, ymm13, ymm1
        vpxor   ymm14, ymm14, ymm2
        vpxor   ymm15, ymm15, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpshufb ymm15, ymm15, ymm8
        vpaddd  ymm8, ymm12, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm13
        vpaddd  ymm10, ymm10, ymm14
        vpaddd  ymm11, ymm11, ymm15
        vpxor   ymm4, ymm4, ymm8
        vpxor   ymm5, ymm5, ymm9
        vpxor   ymm6, ymm6, ymm10
        vpxor   ymm7, ymm7, ymm11
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x1C0]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x40]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x60]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0xE0]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT16+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vmovdqa ymmword ptr [rsp+0x200], ymm8
        vpsrld  ymm8, ymm5, 12
        vpslld  ymm5, ymm5, 20
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 12
        vpslld  ymm6, ymm6, 20
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 12
        vpslld  ymm7, ymm7, 20
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 12
        vpslld  ymm4, ymm4, 20
        vpor    ymm4, ymm4, ymm8
        vpaddd  ymm0, ymm0, ymmword ptr [rsp+0x140]
        vpaddd  ymm1, ymm1, ymmword ptr [rsp+0x180]
        vpaddd  ymm2, ymm2, ymmword ptr [rsp+0x80]
        vpaddd  ymm3, ymm3, ymmword ptr [rsp+0x1A0]
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm1, ymm1, ymm6
        vpaddd  ymm2, ymm2, ymm7
        vpaddd  ymm3, ymm3, ymm4
        vpxor   ymm15, ymm15, ymm0
        vpxor   ymm12, ymm12, ymm1
        vpxor   ymm13, ymm13, ymm2
        vpxor   ymm14, ymm14, ymm3
        vbroadcasti128 ymm8, xmmword ptr [ROT8+rip]
        vpshufb ymm15, ymm15, ymm8
        vpshufb ymm12, ymm12, ymm8
        vpshufb ymm13, ymm13, ymm8
        vpshufb ymm14, ymm14, ymm8
        vpaddd  ymm10, ymm10, ymm15
        vpaddd  ymm11, ymm11, ymm12
        vpaddd  ymm8, ymm13, ymmword ptr [rsp+0x200]
        vpaddd  ymm9, ymm9, ymm14
        vpxor   ymm5, ymm5, ymm10
        vpxor   ymm6, ymm6, ymm11
        vpxor   ymm7, ymm7, ymm8
        vpxor   ymm4, ymm4, ymm9
        vpxor   ymm0, ymm0, ymm8
        vpxor   ymm1, ymm1, ymm9
        vpxor   ymm2, ymm2, ymm10
        vpxor   ymm3, ymm3, ymm11
        vpsrld  ymm8, ymm5, 7
        vpslld  ymm5, ymm5, 25
        vpor    ymm5, ymm5, ymm8
        vpsrld  ymm8, ymm6, 7
        vpslld  ymm6, ymm6, 25
        vpor    ymm6, ymm6, ymm8
        vpsrld  ymm8, ymm7, 7
        vpslld  ymm7, ymm7, 25
        vpor    ymm7, ymm7, ymm8
        vpsrld  ymm8, ymm4, 7
        vpslld  ymm4, ymm4, 25
        vpor    ymm4, ymm4, ymm8
        vpxor   ymm4, ymm4, ymm12
        vpxor   ymm5, ymm5, ymm13
        vpxor   ymm6, ymm6, ymm14
        vpxor   ymm7, ymm7, ymm15
        movzx   eax, byte ptr [rbp+0x38]
        jne     9b
        mov     rbx, qword ptr [rbp+0x50]
        vunpcklps ymm8, ymm0, ymm1
        vunpcklps ymm9, ymm2, ymm3
        vunpckhps ymm10, ymm0, ymm1
        vunpcklps ymm11, ymm4, ymm5
        vunpcklps ymm0, ymm6, ymm7
        vshufps ymm12, ymm8, ymm9, 78
        vblendps ymm1, ymm8, ymm12, 0xCC
        vshufps ymm8, ymm11, ymm0, 78
        vunpckhps ymm13, ymm2, ymm3
        vblendps ymm2, ymm11, ymm8, 0xCC
        vblendps ymm3, ymm12, ymm9, 0xCC
        vperm2f128 ymm12, ymm1, ymm2, 0x20
        vmovups ymmword ptr [rbx], ymm12
        vunpckhps ymm14, ymm4, ymm5
        vblendps ymm4, ymm8, ymm0, 0xCC
        vunpckhps ymm15, ymm6, ymm7
        vperm2f128 ymm7, ymm3, ymm4, 0x20
        vmovups ymmword ptr [rbx+0x20], ymm7
        vshufps ymm5, ymm10, ymm13, 78
        vblendps ymm6, ymm5, ymm13, 0xCC
        vshufps ymm13, ymm14, ymm15, 78
        vblendps ymm10, ymm10, ymm5, 0xCC
        vblendps ymm14, ymm14, ymm13, 0xCC
        vperm2f128 ymm8, ymm10, ymm14, 0x20
        vmovups ymmword ptr [rbx+0x40], ymm8
        vblendps ymm15, ymm13, ymm15, 0xCC
        vperm2f128 ymm13, ymm6, ymm15, 0x20
        vmovups ymmword ptr [rbx+0x60], ymm13
        vperm2f128 ymm9, ymm1, ymm2, 0x31
        vperm2f128 ymm11, ymm3, ymm4, 0x31
        vmovups ymmword ptr [rbx+0x80], ymm9
        vperm2f128 ymm14, ymm10, ymm14, 0x31
        vperm2f128 ymm15, ymm6, ymm15, 0x31
        vmovups ymmword ptr [rbx+0xA0], ymm11
        vmovups ymmword ptr [rbx+0xC0], ymm14
        vmovups ymmword ptr [rbx+0xE0], ymm15
        vmovdqa ymm0, ymmword ptr [rsp+0x220]
        vpaddd  ymm1, ymm0, ymmword ptr [rsp+0x240]
        vmovdqa ymmword ptr [rsp+0x240], ymm1
        vpxor   ymm0, ymm0, ymmword ptr [CMP_MSB_MASK+rip]
        vpxor   ymm2, ymm1, ymmword ptr [CMP_MSB_MASK+rip]
        vpcmpgtd ymm2, ymm0, ymm2
        vmovdqa ymm0, ymmword ptr [rsp+0x260]
        vpsubd  ymm2, ymm0, ymm2
        vmovdqa ymmword ptr [rsp+0x260], ymm2
        add     rdi, 64
        add     rbx, 256
        mov     qword ptr [rbp+0x50], rbx
        sub     rsi, 8
        cmp     rsi, 8
        jnc     2b
        test    rsi, rsi
        jnz     3f
4:
        vzeroupper
        mov     rsp, rbp
        pop     rbp
        pop     rbx
        pop     r12
        pop     r13
        pop     r14
        pop     r15
        ret
.p2align  5
3:
        mov     rbx, qword ptr [rbp+0x50]
        mov     r15, qword ptr [rsp+0x2A0]
        movzx   r13d, byte ptr [rbp+0x38]
        movzx   r12d, byte ptr [rbp+0x48]
        test    rsi, 0x4
        je      3f
        vbroadcasti128 ymm0, xmmword ptr [rcx]
        vbroadcasti128 ymm1, xmmword ptr [rcx+0x10]
        vmovdqa ymm8, ymm0
        vmovdqa ymm9, ymm1
        vbroadcasti128 ymm12, xmmword ptr [rsp+0x240]
        vbroadcasti128 ymm13, xmmword ptr [rsp+0x260]
        vpunpckldq ymm14, ymm12, ymm13
        vpunpckhdq ymm15, ymm12, ymm13
        vpermq  ymm14, ymm14, 0x50
        vpermq  ymm15, ymm15, 0x50
        vbroadcasti128 ymm12, xmmword ptr [BLAKE3_BLOCK_LEN+rip]
        vpblendd ymm14, ymm14, ymm12, 0x44
        vpblendd ymm15, ymm15, ymm12, 0x44
        vmovdqa ymmword ptr [rsp], ymm14
        vmovdqa ymmword ptr [rsp+0x20], ymm15
        mov     r8, qword ptr [rdi]
        mov     r9, qword ptr [rdi+0x8]
        mov     r10, qword ptr [rdi+0x10]
        mov     r11, qword ptr [rdi+0x18]
        movzx   eax, byte ptr [rbp+0x40]
        or      eax, r13d
        xor     edx, edx
.p2align  5
2:
        mov     r14d, eax
        or      eax, r12d
        add     rdx, 64
        cmp     rdx, r15
        cmovne  eax, r14d
        mov     dword ptr [rsp+0x200], eax
        vmovups ymm2, ymmword ptr [r8+rdx-0x40]
        vinsertf128 ymm2, ymm2, xmmword ptr [r9+rdx-0x40], 0x01
        vmovups ymm3, ymmword ptr [r8+rdx-0x30]
        vinsertf128 ymm3, ymm3, xmmword ptr [r9+rdx-0x30], 0x01
        vshufps ymm4, ymm2, ymm3, 136
        vshufps ymm5, ymm2, ymm3, 221
        vmovups ymm2, ymmword ptr [r8+rdx-0x20]
        vinsertf128 ymm2, ymm2, xmmword ptr [r9+rdx-0x20], 0x01
        vmovups ymm3, ymmword ptr [r8+rdx-0x10]
        vinsertf128 ymm3, ymm3, xmmword ptr [r9+rdx-0x10], 0x01
        vshufps ymm6, ymm2, ymm3, 136
        vshufps ymm7, ymm2, ymm3, 221
        vpshufd ymm6, ymm6, 0x93
        vpshufd ymm7, ymm7, 0x93
        vmovups ymm10, ymmword ptr [r10+rdx-0x40]
        vinsertf128 ymm10, ymm10, xmmword ptr [r11+rdx-0x40], 0x01
        vmovups ymm11, ymmword ptr [r10+rdx-0x30]
        vinsertf128 ymm11, ymm11, xmmword ptr [r11+rdx-0x30], 0x01
        vshufps ymm12, ymm10, ymm11, 136
        vshufps ymm13, ymm10, ymm11, 221
        vmovups ymm10, ymmword ptr [r10+rdx-0x20]
        vinsertf128 ymm10, ymm10, xmmword ptr [r11+rdx-0x20], 0x01
        vmovups ymm11, ymmword ptr [r10+rdx-0x10]
        vinsertf128 ymm11, ymm11, xmmword ptr [r11+rdx-0x10], 0x01
        vshufps ymm14, ymm10, ymm11, 136
        vshufps ymm15, ymm10, ymm11, 221
        vpshufd ymm14, ymm14, 0x93
        vpshufd ymm15, ymm15, 0x93
        prefetcht0 [r8+rdx+0x80]
        prefetcht0 [r9+rdx+0x80]
        prefetcht0 [r10+rdx+0x80]
        prefetcht0 [r11+rdx+0x80]
        vpbroadcastd ymm2, dword ptr [rsp+0x200]
        vmovdqa ymm3, ymmword ptr [rsp]
        vmovdqa ymm11, ymmword ptr [rsp+0x20]
        vpblendd ymm3, ymm3, ymm2, 0x88
        vpblendd ymm11, ymm11, ymm2, 0x88
        vbroadcasti128 ymm2, xmmword ptr [BLAKE3_IV+rip]
        vmovdqa ymm10, ymm2
        mov     al, 7
9:
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm8, ymm8, ymm12
        vmovdqa ymmword ptr [rsp+0x40], ymm4
        nop
        vmovdqa ymmword ptr [rsp+0x60], ymm12
        nop
        vpaddd  ymm0, ymm0, ymm1
        vpaddd  ymm8, ymm8, ymm9
        vpxor   ymm3, ymm3, ymm0
        vpxor   ymm11, ymm11, ymm8
        vbroadcasti128 ymm4, xmmword ptr [ROT16+rip]
        vpshufb ymm3, ymm3, ymm4
        vpshufb ymm11, ymm11, ymm4
        vpaddd  ymm2, ymm2, ymm3
        vpaddd  ymm10, ymm10, ymm11
        vpxor   ymm1, ymm1, ymm2
        vpxor   ymm9, ymm9, ymm10
        vpsrld  ymm4, ymm1, 12
        vpslld  ymm1, ymm1, 20
        vpor    ymm1, ymm1, ymm4
        vpsrld  ymm4, ymm9, 12
        vpslld  ymm9, ymm9, 20
        vpor    ymm9, ymm9, ymm4
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm8, ymm8, ymm13
        vpaddd  ymm0, ymm0, ymm1
        vpaddd  ymm8, ymm8, ymm9
        vmovdqa ymmword ptr [rsp+0x80], ymm5
        vmovdqa ymmword ptr [rsp+0xA0], ymm13
        vpxor   ymm3, ymm3, ymm0
        vpxor   ymm11, ymm11, ymm8
        vbroadcasti128 ymm4, xmmword ptr [ROT8+rip]
        vpshufb ymm3, ymm3, ymm4
        vpshufb ymm11, ymm11, ymm4
        vpaddd  ymm2, ymm2, ymm3
        vpaddd  ymm10, ymm10, ymm11
        vpxor   ymm1, ymm1, ymm2
        vpxor   ymm9, ymm9, ymm10
        vpsrld  ymm4, ymm1, 7
        vpslld  ymm1, ymm1, 25
        vpor    ymm1, ymm1, ymm4
        vpsrld  ymm4, ymm9, 7
        vpslld  ymm9, ymm9, 25
        vpor    ymm9, ymm9, ymm4
        vpshufd ymm0, ymm0, 0x93
        vpshufd ymm8, ymm8, 0x93
        vpshufd ymm3, ymm3, 0x4E
        vpshufd ymm11, ymm11, 0x4E
        vpshufd ymm2, ymm2, 0x39
        vpshufd ymm10, ymm10, 0x39
        vpaddd  ymm0, ymm0, ymm6
        vpaddd  ymm8, ymm8, ymm14
        vpaddd  ymm0, ymm0, ymm1
        vpaddd  ymm8, ymm8, ymm9
        vpxor   ymm3, ymm3, ymm0
        vpxor   ymm11, ymm11, ymm8
        vbroadcasti128 ymm4, xmmword ptr [ROT16+rip]
        vpshufb ymm3, ymm3, ymm4
        vpshufb ymm11, ymm11, ymm4
        vpaddd  ymm2, ymm2, ymm3
        vpaddd  ymm10, ymm10, ymm11
        vpxor   ymm1, ymm1, ymm2
        vpxor   ymm9, ymm9, ymm10
        vpsrld  ymm4, ymm1, 12
        vpslld  ymm1, ymm1, 20
        vpor    ymm1, ymm1, ymm4
        vpsrld  ymm4, ymm9, 12
        vpslld  ymm9, ymm9, 20
        vpor    ymm9, ymm9, ymm4
        vpaddd  ymm0, ymm0, ymm7
        vpaddd  ymm8, ymm8, ymm15
        vpaddd  ymm0, ymm0, ymm1
        vpaddd  ymm8, ymm8, ymm9
        vpxor   ymm3, ymm3, ymm0
        vpxor   ymm11, ymm11, ymm8
        vbroadcasti128 ymm4, xmmword ptr [ROT8+rip]
        vpshufb ymm3, ymm3, ymm4
        vpshufb ymm11, ymm11, ymm4
        vpaddd  ymm2, ymm2, ymm3
        vpaddd  ymm10, ymm10, ymm11
        vpxor   ymm1, ymm1, ymm2
        vpxor   ymm9, ymm9, ymm10
        vpsrld  ymm4, ymm1, 7
        vpslld  ymm1, ymm1, 25
        vpor    ymm1, ymm1, ymm4
        vpsrld  ymm4, ymm9, 7
        vpslld  ymm9, ymm9, 25
        vpor    ymm9, ymm9, ymm4
        vpshufd ymm0, ymm0, 0x39
        vpshufd ymm8, ymm8, 0x39
        vpshufd ymm3, ymm3, 0x4E
        vpshufd ymm11, ymm11, 0x4E
        vpshufd ymm2, ymm2, 0x93
        vpshufd ymm10, ymm10, 0x93
        dec     al
        je      9f
        vmovdqa ymm4, ymmword ptr [rsp+0x40]
        vmovdqa ymm5, ymmword ptr [rsp+0x80]
        vshufps ymm12, ymm4, ymm5, 214
        vpshufd ymm13, ymm4, 0x0F
        vpshufd ymm4, ymm12, 0x39
        vshufps ymm12, ymm6, ymm7, 250
        vpblendd ymm13, ymm13, ymm12, 0xAA
        vpunpcklqdq ymm12, ymm7, ymm5
        vpblendd ymm12, ymm12, ymm6, 0x88
        vpshufd ymm12, ymm12, 0x78
        vpunpckhdq ymm5, ymm5, ymm7
        vpunpckldq ymm6, ymm6, ymm5
        vpshufd ymm7, ymm6, 0x1E
        vmovdqa ymmword ptr [rsp+0x40], ymm13
        vmovdqa ymmword ptr [rsp+0x80], ymm12
        vmovdqa ymm12, ymmword ptr [rsp+0x60]
        vmovdqa ymm13, ymmword ptr [rsp+0xA0]
        vshufps ymm5, ymm12, ymm13, 214
        vpshufd ymm6, ymm12, 0x0F
        vpshufd ymm12, ymm5, 0x39
        vshufps ymm5, ymm14, ymm15, 250
        vpblendd ymm6, ymm6, ymm5, 0xAA
        vpunpcklqdq ymm5, ymm15, ymm13
        vpblendd ymm5, ymm5, ymm14, 0x88
        vpshufd ymm5, ymm5, 0x78
        vpunpckhdq ymm13, ymm13, ymm15
        vpunpckldq ymm14, ymm14, ymm13
        vpshufd ymm15, ymm14, 0x1E
        vmovdqa ymm13, ymm6
        vmovdqa ymm14, ymm5
        vmovdqa ymm5, ymmword ptr [rsp+0x40]
        vmovdqa ymm6, ymmword ptr [rsp+0x80]
        jmp     9b
9:
        vpxor   ymm0, ymm0, ymm2
        vpxor   ymm1, ymm1, ymm3
        vpxor   ymm8, ymm8, ymm10
        vpxor   ymm9, ymm9, ymm11
        mov     eax, r13d
        cmp     rdx, r15
        jne     2b
        vmovdqu xmmword ptr [rbx], xmm0
        vmovdqu xmmword ptr [rbx+0x10], xmm1
        vextracti128 xmmword ptr [rbx+0x20], ymm0, 0x01
        vextracti128 xmmword ptr [rbx+0x30], ymm1, 0x01
        vmovdqu xmmword ptr [rbx+0x40], xmm8
        vmovdqu xmmword ptr [rbx+0x50], xmm9
        vextracti128 xmmword ptr [rbx+0x60], ymm8, 0x01
        vextracti128 xmmword ptr [rbx+0x70], ymm9, 0x01
        vmovaps xmm8, xmmword ptr [rsp+0x280]
        vmovaps xmm0, xmmword ptr [rsp+0x240]
        vmovaps xmm1, xmmword ptr [rsp+0x250]
        vmovaps xmm2, xmmword ptr [rsp+0x260]
        vmovaps xmm3, xmmword ptr [rsp+0x270]
        vblendvps xmm0, xmm0, xmm1, xmm8
        vblendvps xmm2, xmm2, xmm3, xmm8
        vmovaps xmmword ptr [rsp+0x240], xmm0
        vmovaps xmmword ptr [rsp+0x260], xmm2
        add     rbx, 128
        add     rdi, 32
        sub     rsi, 4
3:
        test    rsi, 0x2
        je      3f
        vbroadcasti128 ymm0, xmmword ptr [rcx]
        vbroadcasti128 ymm1, xmmword ptr [rcx+0x10]
        vmovd   xmm13, dword ptr [rsp+0x240]
        vpinsrd xmm13, xmm13, dword ptr [rsp+0x260], 1
        vpinsrd xmm13, xmm13, dword ptr [BLAKE3_BLOCK_LEN+rip], 2
        vmovd   xmm14, dword ptr [rsp+0x244]
        vpinsrd xmm14, xmm14, dword ptr [rsp+0x264], 1
        vpinsrd xmm14, xmm14, dword ptr [BLAKE3_BLOCK_LEN+rip], 2
        vinserti128 ymm13, ymm13, xmm14, 0x01
        vbroadcasti128 ymm14, xmmword ptr [ROT16+rip]
        vbroadcasti128 ymm15, xmmword ptr [ROT8+rip]
        mov     r8, qword ptr [rdi]
        mov     r9, qword ptr [rdi+0x8]
        movzx   eax, byte ptr [rbp+0x40]
        or      eax, r13d
        xor     edx, edx
.p2align  5
2:
        mov     r14d, eax
        or      eax, r12d
        add     rdx, 64
        cmp     rdx, r15
        cmovne  eax, r14d
        mov     dword ptr [rsp+0x200], eax
        vbroadcasti128 ymm2, xmmword ptr [BLAKE3_IV+rip]
        vpbroadcastd ymm8, dword ptr [rsp+0x200]
        vpblendd ymm3, ymm13, ymm8, 0x88
        vmovups ymm8, ymmword ptr [r8+rdx-0x40]
        vinsertf128 ymm8, ymm8, xmmword ptr [r9+rdx-0x40], 0x01
        vmovups ymm9, ymmword ptr [r8+rdx-0x30]
        vinsertf128 ymm9, ymm9, xmmword ptr [r9+rdx-0x30], 0x01
        vshufps ymm4, ymm8, ymm9, 136
        vshufps ymm5, ymm8, ymm9, 221
        vmovups ymm8, ymmword ptr [r8+rdx-0x20]
        vinsertf128 ymm8, ymm8, xmmword ptr [r9+rdx-0x20], 0x01
        vmovups ymm9, ymmword ptr [r8+rdx-0x10]
        vinsertf128 ymm9, ymm9, xmmword ptr [r9+rdx-0x10], 0x01
        vshufps ymm6, ymm8, ymm9, 136
        vshufps ymm7, ymm8, ymm9, 221
        vpshufd ymm6, ymm6, 0x93
        vpshufd ymm7, ymm7, 0x93
        mov     al, 7
9:
        vpaddd  ymm0, ymm0, ymm4
        vpaddd  ymm0, ymm0, ymm1
        vpxor   ymm3, ymm3, ymm0
        vpshufb ymm3, ymm3, ymm14
        vpaddd  ymm2, ymm2, ymm3
        vpxor   ymm1, ymm1, ymm2
        vpsrld  ymm8, ymm1, 12
        vpslld  ymm1, ymm1, 20
        vpor    ymm1, ymm1, ymm8
        vpaddd  ymm0, ymm0, ymm5
        vpaddd  ymm0, ymm0, ymm1
        vpxor   ymm3, ymm3, ymm0
        vpshufb ymm3, ymm3, ymm15
        vpaddd  ymm2, ymm2, ymm3
        vpxor   ymm1, ymm1, ymm2
        vpsrld  ymm8, ymm1, 7
        vpslld  ymm1, ymm1, 25
        vpor    ymm1, ymm1, ymm8
        vpshufd ymm0, ymm0, 0x93
        vpshufd ymm3, ymm3, 0x4E
        vpshufd ymm2, ymm2, 0x39
        vpaddd  ymm0, ymm0, ymm6
        vpaddd  ymm0, ymm0, ymm1
        vpxor   ymm3, ymm3, ymm0
        vpshufb ymm3, ymm3, ymm14
        vpaddd  ymm2, ymm2, ymm3
        vpxor   ymm1, ymm1, ymm2
        vpsrld  ymm8, ymm1, 12
        vpslld  ymm1, ymm1, 20
        vpor    ymm1, ymm1, ymm8
        vpaddd  ymm0, ymm0, ymm7
        vpaddd  ymm0, ymm0, ymm1
        vpxor   ymm3, ymm3, ymm0
        vpshufb ymm3, ymm3, ymm15
        vpaddd  ymm2, ymm2, ymm3
        vpxor   ymm1, ymm1, ymm2
        vpsrld  ymm8, ymm1, 7
        vpslld  ymm1, ymm1, 25
        vpor    ymm1, ymm1, ymm8
        vpshufd ymm0, ymm0, 0x39
        vpshufd ymm3, ymm3, 0x4E
        vpshufd ymm2, ymm2, 0x93
        dec     al
        jz      9f
        vshufps ymm8, ymm4, ymm5, 214
        vpshufd ymm9, ymm4, 0x0F
        vpshufd ymm4, ymm8, 0x39
        vshufps ymm8, ymm6, ymm7, 250
        vpblendd ymm9, ymm9, ymm8, 0xAA
        vpunpcklqdq ymm8, ymm7, ymm5
        vpblendd ymm8, ymm8, ymm6, 0x88
        vpshufd ymm8, ymm8, 0x78
        vpunpckhdq ymm5, ymm5, ymm7
        vpunpckldq ymm6, ymm6, ymm5
        vpshufd ymm7, ymm6, 0x1E
        vmovdqa ymm5, ymm9
        vmovdqa ymm6, ymm8
        jmp     9b
9:
        vpxor   ymm0, ymm0, ymm2
        vpxor   ymm1, ymm1, ymm3
        mov     eax, r13d
        cmp     rdx, r15
        jne     2b
        vmovdqu xmmword ptr [rbx], xmm0
        vmovdqu xmmword ptr [rbx+0x10], xmm1
        vextracti128 xmmword ptr [rbx+0x20], ymm0, 0x01
        vextracti128 xmmword ptr [rbx+0x30], ymm1, 0x01
        vmovaps ymm8, ymmword ptr [rsp+0x280]
        vmovaps ymm0, ymmword ptr [rsp+0x240]
        vmovups ymm1, ymmword ptr [rsp+0x248]
        vmovaps ymm2, ymmword ptr [rsp+0x260]
        vmovups ymm3, ymmword ptr [rsp+0x268]
        vblendvps ymm0, ymm0, ymm1, ymm8
        vblendvps ymm2, ymm2, ymm3, ymm8
        vmovaps ymmword ptr [rsp+0x240], ymm0
        vmovaps ymmword ptr [rsp+0x260], ymm2
        add     rbx, 64
        add     rdi, 16
        sub     rsi, 2
3:
        test    rsi, 0x1
        je      4b
        vmovdqu xmm0, xmmword ptr [rcx]
        vmovdqu xmm1, xmmword ptr [rcx+0x10]
        vmovd   xmm3, dword ptr [rsp+0x240]
        vpinsrd xmm3, xmm3, dword ptr [rsp+0x260], 1
        vpinsrd xmm13, xmm3, dword ptr [BLAKE3_BLOCK_LEN+rip], 2
        vmovdqa xmm14, xmmword ptr [ROT16+rip]
        vmovdqa xmm15, xmmword ptr [ROT8+rip]
        mov     r8, qword ptr [rdi]
        movzx   eax, byte ptr [rbp+0x40]
        or      eax, r13d
        xor     edx, edx
.p2align  5
2:
        mov     r14d, eax
        or      eax, r12d
        add     rdx, 64
        cmp     rdx, r15
        cmovne  eax, r14d
        vmovdqa xmm2, xmmword ptr [BLAKE3_IV+rip]
        vmovdqa xmm3, xmm13
        vpinsrd xmm3, xmm3, eax, 3
        vmovups xmm8, xmmword ptr [r8+rdx-0x40]
        vmovups xmm9, xmmword ptr [r8+rdx-0x30]
        vshufps xmm4, xmm8, xmm9, 136
        vshufps xmm5, xmm8, xmm9, 221
        vmovups xmm8, xmmword ptr [r8+rdx-0x20]
        vmovups xmm9, xmmword ptr [r8+rdx-0x10]
        vshufps xmm6, xmm8, xmm9, 136
        vshufps xmm7, xmm8, xmm9, 221
        vpshufd xmm6, xmm6, 0x93
        vpshufd xmm7, xmm7, 0x93
        mov     al, 7
9:
        vpaddd  xmm0, xmm0, xmm4
        vpaddd  xmm0, xmm0, xmm1
        vpxor   xmm3, xmm3, xmm0
        vpshufb xmm3, xmm3, xmm14
        vpaddd  xmm2, xmm2, xmm3
        vpxor   xmm1, xmm1, xmm2
        vpsrld  xmm8, xmm1, 12
        vpslld  xmm1, xmm1, 20
        vpor    xmm1, xmm1, xmm8
        vpaddd  xmm0, xmm0, xmm5
        vpaddd  xmm0, xmm0, xmm1
        vpxor   xmm3, xmm3, xmm0
        vpshufb xmm3, xmm3, xmm15
        vpaddd  xmm2, xmm2, xmm3
        vpxor   xmm1, xmm1, xmm2
        vpsrld  xmm8, xmm1, 7
        vpslld  xmm1, xmm1, 25
        vpor    xmm1, xmm1, xmm8
        vpshufd xmm0, xmm0, 0x93
        vpshufd xmm3, xmm3, 0x4E
        vpshufd xmm2, xmm2, 0x39
        vpaddd  xmm0, xmm0, xmm6
        vpaddd  xmm0, xmm0, xmm1
        vpxor   xmm3, xmm3, xmm0
        vpshufb xmm3, xmm3, xmm14
        vpaddd  xmm2, xmm2, xmm3
        vpxor   xmm1, xmm1, xmm2
        vpsrld  xmm8, xmm1, 12
        vpslld  xmm1, xmm1, 20
        vpor    xmm1, xmm1, xmm8
        vpaddd  xmm0, xmm0, xmm7
        vpaddd  xmm0, xmm0, xmm1
        vpxor   xmm3, xmm3, xmm0
        vpshufb xmm3, xmm3, xmm15
        vpaddd  xmm2, xmm2, xmm3
        vpxor   xmm1, xmm1, xmm2
        vpsrld  xmm8, xmm1, 7
        vpslld  xmm1, xmm1, 25
        vpor    xmm1, xmm1, xmm8
        vpshufd xmm0, xmm0, 0x39
        vpshufd xmm3, xmm3, 0x4E
        vpshufd xmm2, xmm2, 0x93
        dec     al
        jz      9f
        vshufps xmm8, xmm4, xmm5, 214
        vpshufd xmm9, xmm4, 0x0F
        vpshufd xmm4, xmm8, 0x39
        vshufps xmm8, xmm6, xmm7, 250
        vpblendd xmm9, xmm9, xmm8, 0xAA
        vpunpcklqdq xmm8, xmm7, xmm5
        vpblendd xmm8, xmm8, xmm6, 0x88
        vpshufd xmm8, xmm8, 0x78
        vpunpckhdq xmm5, xmm5, xmm7
        vpunpckldq xmm6, xmm6, xmm5
        vpshufd xmm7, xmm6, 0x1E
        vmovdqa xmm5, xmm9
        vmovdqa xmm6, xmm8
        jmp     9b
9:
        vpxor   xmm0, xmm0, xmm2
        vpxor   xmm1, xmm1, xmm3
        mov     eax, r13d
        cmp     rdx, r15
        jne     2b
        vmovdqu xmmword ptr [rbx], xmm0
        vmovdqu xmmword ptr [rbx+0x10], xmm1
        jmp     4b


.extern blake3_compress_xof
.global _blake3_xof_many_avx22
.global blake3_xof_many_avx22

# ymmX: stateX
# RDI: cv
# RDX: block_len
# RCX: counter/counter_low
# RBX: counter_high
# R8: flags
.macro init_state
        vmovdqa ymm15, ymmword ptr [rip+ADD0] # ymm15: inc
        vmovd xmm12, rcx
        vpbroadcastd ymm12, xmm12             # ymm12: counter_low
        mov rbx, rcx                          # rbx: counter_high
        vpbroadcastd ymm0, dword ptr [rdi]    # ymmX: stateX (from cv[X])
        vpbroadcastd ymm1, dword ptr [rdi+4]  # ymmX: stateX (from cv[X])
        vpbroadcastd ymm2, dword ptr [rdi+8]  # ymmX: stateX (from cv[X])
        vpbroadcastd ymm3, dword ptr [rdi+12] # ymmX: stateX (from cv[X])
        vpbroadcastd ymm4, dword ptr [rdi+16] # ymmX: stateX (from cv[X])
        shr rbx, 32
        vpbroadcastd ymm5, dword ptr [rdi+20] # ymmX: stateX (from cv[X])
        vpbroadcastd ymm6, dword ptr [rdi+24] # ymmX: stateX (from cv[X])
        vpbroadcastd ymm7, dword ptr [rdi+28] # ymmX: stateX (from cv[X])
        vpbroadcastd ymm8, dword ptr [rip+BLAKE3_IV]      # ymmX: stateX (from IV[X-8])
        vpbroadcastd ymm9, dword ptr [rip+BLAKE3_IV+4]    # ymmX: stateX (from IV[X-8])
        vpbroadcastd ymm10, dword ptr [rip+BLAKE3_IV+8]   # ymmX: stateX (from IV[X-8])
        vpbroadcastd ymm11, dword ptr [rip+BLAKE3_IV+12]  # ymmX: stateX (from IV[X-8])
        # Only works if the upper part doesn't overflow
        vpaddd ymm12, ymm12, ymm15                        # ymm12: inc by 0, 1, ..., 7
        vmovd xmm13, rbx
        vmovd xmm14, rdx
        vmovd xmm15, r8
        vpbroadcastd ymm13, xmm13                         # ymm13: counter_high
        vpbroadcastd ymm14, xmm14                         # ymm14: block_len
        vpbroadcastd ymm15, xmm15                         # ymm15: flags
.endm

# Inputs: v0, v2. Temporaries: v1, v3
.macro rot2 v0, v1, imm1, \
            v2, v3, imm2
        .set _imm1, \imm1
        .set _imm2, \imm2
        vpsrld \v1, \v0, _imm1
        vpsrld \v3, \v2, _imm2
        vpslld \v0, \v0, 32-_imm1
        vpslld \v2, \v2, 32-_imm2
        vpor \v0, \v0, \v1
        vpor \v2, \v2, \v3
.endm

.macro g2_avx2 vA, vB, vC, vD, vx, vy, \
               vA2, vB2, vC2, vD2, vx2, vy2
        vpaddd \vA, \vA, \vB
        vpaddd \vA2, \vA2, \vB2
        vpaddd \vA, \vA, \vx
        vpaddd \vA2, \vA2, \vx2
        vpxor \vD, \vA, \vD
        vpxor \vD2, \vA2, \vD2
        rot2 \vD, \vx, 16, \
             \vD2, \vx2, 16

        vpaddd \vC, \vC, \vD
        vpaddd \vC2, \vC2, \vD2
        vpxor \vB, \vC, \vB
        vpxor \vB2, \vC2, \vB2
        rot2 \vB, \vx, 12, \
             \vB2, \vx2, 12

        vpaddd \vA, \vA, \vy
        vpaddd \vA2, \vA2, \vy2
        vpaddd \vA, \vA, \vB
        vpaddd \vA2, \vA2, \vB2
        vpxor \vD, \vA, \vD
        vpxor \vD2, \vA2, \vD2
        rot2 \vD, \vy, 8, \
             \vD2, \vy2, 8

        vpaddd \vC, \vC, \vD
        vpaddd \vC2, \vC2, \vD2
        vpxor \vB, \vC, \vB
        vpxor \vB2, \vC2, \vB2
        rot2 \vB, \vy, 7, \
             \vB2, \vy2, 7
.endm

# ymmX: stateX
# RSI: msg
# r12-15: temporaries
# 256-byte allocated in stack
.macro round_fn round
        .set _round, \round
        # Mix columns 0 and 1
        vmovdqu [rsp], ymm10                            # save 4 vectors on the stack
        vmovdqu [rsp+32], ymm11
        vmovdqu [rsp+64], ymm14
        vmovdqu [rsp+96], ymm15
        movzx r12, byte ptr [rip+SCHEDULE+16*_round]    # Load msg[schedule[0-3]]
        movzx r13, byte ptr [rip+SCHEDULE+16*_round+1]
        movzx r14, byte ptr [rip+SCHEDULE+16*_round+2]
        movzx r15, byte ptr [rip+SCHEDULE+16*_round+3]
        vpbroadcastd ymm10, dword ptr [rsi + r12*4]
        vpbroadcastd ymm11, dword ptr [rsi + r13*4]
        vpbroadcastd ymm14, dword ptr [rsi + r14*4]
        vpbroadcastd ymm15, dword ptr [rsi + r15*4]
        g2_avx2 ymm0, ymm4, ymm8, ymm12, ymm10, ymm11, \
                ymm1, ymm5, ymm9, ymm13, ymm14, ymm15

        # Mix columns 2 and 3
        vmovdqu [rsp+128], ymm8                         # save 4 vectors on the stack
        vmovdqu [rsp+160], ymm13
        vmovdqu [rsp+192], ymm4
        vmovdqu [rsp+224], ymm9
        vmovdqu ymm10, [rsp]                            # restore 4 vectors
        vmovdqu ymm11, [rsp+32]
        vmovdqu ymm14, [rsp+64]
        vmovdqu ymm15, [rsp+96]
        movzx r12, byte ptr [rip+SCHEDULE+16*_round+4]  # Load msg[schedule[4-7]]
        movzx r13, byte ptr [rip+SCHEDULE+16*_round+5]
        movzx r14, byte ptr [rip+SCHEDULE+16*_round+6]
        movzx r15, byte ptr [rip+SCHEDULE+16*_round+7]
        vpbroadcastd ymm8, dword ptr [rsi + r12*4]
        vpbroadcastd ymm13, dword ptr [rsi + r13*4]
        vpbroadcastd ymm4, dword ptr [rsi + r14*4]
        vpbroadcastd ymm9, dword ptr [rsi + r15*4]
        g2_avx2 ymm2, ymm6, ymm10, ymm14, ymm8, ymm13, \
                ymm3, ymm7, ymm11, ymm15, ymm4, ymm9

        # Mix diagonals 0 and 1
        movzx r12, byte ptr [rip+SCHEDULE+16*_round+8]  # Load msg[schedule[8-11]]
        movzx r13, byte ptr [rip+SCHEDULE+16*_round+9]
        movzx r14, byte ptr [rip+SCHEDULE+16*_round+10]
        movzx r15, byte ptr [rip+SCHEDULE+16*_round+11]
        vpbroadcastd ymm8, dword ptr [rsi + r12*4]
        vpbroadcastd ymm13, dword ptr [rsi + r13*4]
        vpbroadcastd ymm4, dword ptr [rsi + r14*4]
        vpbroadcastd ymm9, dword ptr [rsi + r15*4]
        g2_avx2 ymm0, ymm5, ymm10, ymm15, ymm8, ymm13, \
                ymm1, ymm6, ymm11, ymm12, ymm4, ymm9

        # Mix diagonals 2 and 3
        vmovdqu [rsp], ymm0                             # save 4 vectors on the stack
        vmovdqu [rsp+32], ymm1
        vmovdqu [rsp+64], ymm5
        vmovdqu [rsp+96], ymm6
        vmovdqu ymm8, [rsp+128]                         # restore 4 vectors
        vmovdqu ymm13, [rsp+160]
        vmovdqu ymm4, [rsp+192]
        vmovdqu ymm9, [rsp+224]
        movzx r12, byte ptr [rip+SCHEDULE+16*_round+12] # Load msg[schedule[12-15]]
        movzx r13, byte ptr [rip+SCHEDULE+16*_round+13]
        movzx r14, byte ptr [rip+SCHEDULE+16*_round+14]
        movzx r15, byte ptr [rip+SCHEDULE+16*_round+15]
        vpbroadcastd ymm8, dword ptr [rsi + r12*4]
        vpbroadcastd ymm13, dword ptr [rsi + r13*4]
        vpbroadcastd ymm4, dword ptr [rsi + r14*4]
        vpbroadcastd ymm9, dword ptr [rsi + r15*4]
        g2_avx2 ymm2, ymm7, ymm8, ymm13, ymm0, ymm1, \
                ymm3, ymm4, ymm9, ymm14, ymm5, ymm6

        vmovdqu ymm0, [rsp]                             # restore 4 vectors
        vmovdqu ymm1, [rsp+32]
        vmovdqu ymm5, [rsp+64]
        vmovdqu ymm6, [rsp+96]
.endm

# v0-v7: stateX
# v15: tmp
# R9: out
.macro transpose_store v0, v1, v2, v3, \
                       v4, v5, v6, v7, \
                       v15, offset
        .set _off, \offset
        # Tranpose Logic
        vmovdqu [rsp], \v15                            # save 1 vector on the stack
        # First Layer
        vmovdqa \v15, \v0
        vpunpckldq \v0, \v15, \v1
        vpunpckhdq \v1, \v15, \v1
        vmovdqa \v15, \v2
        vpunpckldq \v2, \v15, \v3
        vpunpckhdq \v3, \v15, \v3
        vmovdqa \v15, \v4
        vpunpckldq \v4, \v15, \v5
        vpunpckhdq \v5, \v15, \v5
        vmovdqa \v15, \v6
        vpunpckldq \v6, \v15, \v7
        vpunpckhdq \v7, \v15, \v7
        # Second Layer
        vmovdqa \v15, \v0
        vpunpckldq \v0, \v15, \v2
        vpunpckhdq \v2, \v15, \v2
        vmovdqa \v15, \v1
        vpunpckldq \v1, \v15, \v3
        vpunpckhdq \v3, \v15, \v3
        vmovdqa \v15, \v4
        vpunpckldq \v4, \v15, \v6
        vpunpckhdq \v6, \v15, \v6
        vmovdqa \v15, \v5
        vpunpckldq \v5, \v15, \v7
        vpunpckhdq \v7, \v15, \v7
        # Third Layer
        vmovdqa \v15, \v0
        vperm2i128 \v0, \v15, \v4, 0x20
        vperm2i128 \v4, \v15, \v4, 0x31
        vmovdqa \v15, \v1
        vperm2i128 \v1, \v15, \v5, 0x20
        vperm2i128 \v5, \v15, \v5, 0x31
        vmovdqa \v15, \v2
        vperm2i128 \v2, \v15, \v6, 0x20
        vperm2i128 \v6, \v15, \v6, 0x31
        vmovdqa \v15, \v3
        vperm2i128 \v3, \v15, \v7, 0x20
        vperm2i128 \v7, \v15, \v7, 0x31
        # Store
        vmovdqu [r9+_off], \v0
        vmovdqu [r9+_off+64], \v2
        vmovdqu [r9+_off+128], \v1
        vmovdqu [r9+_off+192], \v3
        vmovdqu [r9+_off+256], \v4
        vmovdqu [r9+_off+320], \v6
        vmovdqu [r9+_off+384], \v5
        vmovdqu [r9+_off+448], \v7
.endm

# ymmX: stateX
# R9: out
.macro finalize_state_0
        # XOR Logic
        vpxor ymm0, ymm0, ymm8
        vpxor ymm1, ymm1, ymm9
        vpxor ymm2, ymm2, ymm10
        vpxor ymm3, ymm3, ymm11
        vpxor ymm4, ymm4, ymm12
        vpxor ymm5, ymm5, ymm13
        vpxor ymm6, ymm6, ymm14
        vpxor ymm7, ymm7, ymm15
        # Tranpose Logic
        transpose_store ymm0, ymm1, ymm2, ymm3, \
                        ymm4, ymm5, ymm6, ymm7, \
                        ymm15, 0
.endm

# ymmX: stateX
# RDI: cv
# R9: out
.macro finalize_state_1
        vmovdqu ymm15, [rsp]                  # restore 1 vector on the stack
        # XOR Logic
        vpbroadcastd ymm0, dword ptr [rdi]    # ymmX: stateX (from cv[X])
        vpbroadcastd ymm1, dword ptr [rdi+4]  # ymmX: stateX (from cv[X])
        vpbroadcastd ymm2, dword ptr [rdi+8]  # ymmX: stateX (from cv[X])
        vpbroadcastd ymm3, dword ptr [rdi+12] # ymmX: stateX (from cv[X])
        vpbroadcastd ymm4, dword ptr [rdi+16] # ymmX: stateX (from cv[X])
        vpbroadcastd ymm5, dword ptr [rdi+20] # ymmX: stateX (from cv[X])
        vpbroadcastd ymm6, dword ptr [rdi+24] # ymmX: stateX (from cv[X])
        vpbroadcastd ymm7, dword ptr [rdi+28] # ymmX: stateX (from cv[X])
        vpxor ymm8, ymm8, ymm0
        vpxor ymm9, ymm9, ymm1
        vpxor ymm10, ymm10, ymm2
        vpxor ymm11, ymm11, ymm3
        vpxor ymm12, ymm12, ymm4
        vpxor ymm13, ymm13, ymm5
        vpxor ymm14, ymm14, ymm6
        vpxor ymm15, ymm15, ymm7
        # Tranpose Logic
        transpose_store ymm8, ymm9, ymm10, ymm11, \
                        ymm12, ymm13, ymm14, ymm15, \
                        ymm0, 32
.endm

# RDI: cv
# RSI: block
# RDX: block_len
# RCX: counter
# R8: flags
# R9: out
# Stack (RSP) = outblocks
_blake3_xof_many_avx22:
blake3_xof_many_avx22:
        mov r10, [rsp + 8]      # outblocks (r10)
        push rbx                # save rbx, r12-15 on the stack
        push r12
        push r13
        push r14
        push r15
        mov r11, r10            # outblocks8 (r10): r11 -> temp
        and r11, 7
        sub r10, r11
        test r10, r10
        jz remaining_xof
        mov r11, rcx            # iterate until counter + outblocks8 (r11)
        add r11, r10
        sub rsp, 256            # reserve stack space for 8 vector register
xof_8_blocks:
        # Init state
        init_state
        # 7 Rounds
        round_fn 0
        round_fn 1
        round_fn 2
        round_fn 3
        round_fn 4
        round_fn 5
        round_fn 6
        # Finalize State
        finalize_state_0
        finalize_state_1
        # Loop Logic
        add rcx, 8
        add r9, 64
        cmp rcx, r11
        jne  xof_8_blocks       # jump if counter < counter + outblocks8 -> Restart the loop
        add rsp, 256
remaining_xof:
        mov rbx, [rsp + 48]     # outblocks (rbx)
        sub rbx, r10            # i = outblocks - outblocks8
        test rbx, rbx
        jz blake3_avx2_xof_end
        mov rbx, rdi
        mov r13, rsi
        mov r14, rdx
        push rcx
        mov r15, r8
        push r9
xof_1_block:
        mov rdi, r12            # caller-saved
        mov rsi, r13
        mov rdx, r14
        mov r8, r15
        call    blake3_compress_xof
        pop rcx                 # caller-saved + inc logic
        pop r9
        add rcx, 1
        add r9, 64
        sub rbx, 1
        push rcx
        push r9
        test rbx, rbx           # continue if i > 0
        jnz xof_1_block
blake3_avx2_xof_end:
        pop r15
        pop r14
        pop r13
        pop r12
        pop rbx
        ret

#ifdef __APPLE__
.static_data
#else
.section .rodata
#endif
.p2align  6
ADD0:
        .long  0, 1, 2, 3, 4, 5, 6, 7
ADD1:
        .long  8, 8, 8, 8, 8, 8, 8, 8
BLAKE3_IV_0:
        .long  0x6A09E667, 0x6A09E667, 0x6A09E667, 0x6A09E667
        .long  0x6A09E667, 0x6A09E667, 0x6A09E667, 0x6A09E667
BLAKE3_IV_1:
        .long  0xBB67AE85, 0xBB67AE85, 0xBB67AE85, 0xBB67AE85
        .long  0xBB67AE85, 0xBB67AE85, 0xBB67AE85, 0xBB67AE85
BLAKE3_IV_2:
        .long  0x3C6EF372, 0x3C6EF372, 0x3C6EF372, 0x3C6EF372
        .long  0x3C6EF372, 0x3C6EF372, 0x3C6EF372, 0x3C6EF372
BLAKE3_IV_3:
        .long  0xA54FF53A, 0xA54FF53A, 0xA54FF53A, 0xA54FF53A
        .long  0xA54FF53A, 0xA54FF53A, 0xA54FF53A, 0xA54FF53A
BLAKE3_BLOCK_LEN:
        .long  0x00000040, 0x00000040, 0x00000040, 0x00000040
        .long  0x00000040, 0x00000040, 0x00000040, 0x00000040
ROT16:
        .byte  2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13
ROT8:
        .byte  1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12
CMP_MSB_MASK:
        .long  0x80000000, 0x80000000, 0x80000000, 0x80000000
        .long  0x80000000, 0x80000000, 0x80000000, 0x80000000
BLAKE3_IV:
        .long  0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A
        .long  0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19
SCHEDULE:
        .byte  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
        .byte  2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8
        .byte  3, 4, 10, 12, 13, 2, 7, 14, 6, 5, 9, 0, 11, 15, 8, 1
        .byte  10, 7, 12, 9, 14, 3, 13, 15, 4, 0, 11, 2, 5, 8, 1, 6
        .byte  12, 13, 9, 11, 15, 10, 14, 8, 7, 2, 5, 3, 0, 1, 6, 4
        .byte  9, 14, 11, 5, 8, 12, 15, 1, 13, 3, 0, 10, 2, 6, 4, 7
        .byte  11, 15, 5, 0, 1, 9, 8, 6, 14, 10, 2, 12, 3, 4, 7, 13
